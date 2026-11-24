// ==============================================================================
// [STAGE 7b] Validation of the atomic relaxation after photoabsorption in lead.
//
// Data: ENDF/B-VIII.0 (EPICS2017) for lead in data/: atom-082_Pb_000.endf
// (relaxation), photoat-082_Pb_000.endf (subshell photoelectric cross sections).
//
//   TEST 1  Shell data: binding energies vs X-ray Data Booklet, K fluorescence
//           yield vs Krause, K line energies vs Bearden (measured)
//   TEST 2  Subshell selection: ENDF shell fractions at K and L edges vs NIST
//           XCOM edge jumps, and sampled vs computed frequencies
//   TEST 3  Cascade sampling: K X-rays per K vacancy (= omega_K) and line fractions
//   TEST 4  Energy conservation of photoabsorption + relaxation, event by event
//   TEST 5  K X-ray escape from a thick lead block vs the analytic formula
//
// Run:  ./ValidationRelaxApp              (~5.4 s)
//       ./ValidationRelaxApp quick        (fewer events)
//       ./ValidationRelaxApp escape 7     (TEST 5 only, seed 7)
// ==============================================================================

#include "AtomicDeexcitation.h"
#include "AtomicTransitionManager.h"
#include "Box.h"
#include "Ionization.h"
#include "Material.h"
#include "Particle.h"
#include "PhotoElectricEffect.h"
#include "PhotoElectricShellData.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

using PhysicalConstants::m_e;

namespace {

const Particle gammaDef(22, "gamma", 0.0, 0.0, 1.0, -1.0);
const Particle electronDef(11, "e-", m_e, -1.0, 0.5, -1.0);
constexpr int kZ = 82;
constexpr double kMaxPull = 3.0;

std::shared_ptr<Material> makeLead() {
    const DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

void printHeader(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

int findShell(const AtomicTransitionManager& atom, const std::string& name) {
    for (int i = 0; i < atom.getNumberOfShells(); ++i) {
        if (AtomicTransitionManager::getShellName(atom.getShell(i).designator) == name) {
            return i;
        }
    }
    throw std::runtime_error("no shell " + name);
}

int findPhotoShell(const PhotoElectricShellData& pe, const std::string& name) {
    for (int i = 0; i < pe.getNumberOfShells(); ++i) {
        if (AtomicTransitionManager::getShellName(pe.getDesignator(i)) == name) {
            return i;
        }
    }
    throw std::runtime_error("no photoelectric shell " + name);
}

// Probability and energy of the radiative transition origin -> vacancy (0 if absent)
AtomicTransitionManager::RadiativeTransition findLine(const AtomicTransitionManager& atom, int vacancy, int origin) {
    for (const auto& t : atom.getShell(vacancy).radiative) {
        if (t.originShell == origin) {
            return t;
        }
    }
    return {origin, 0.0, 0.0, 0.0};
}

// ==============================================================================
// TEST 1  Shell data vs reference values
// ==============================================================================
bool testShellData() {
    printHeader("TEST 1  Lead shell data (ENDF/B-VIII.0) vs reference values");
    const AtomicTransitionManager& eadl = AtomicTransitionManager::get(kZ);
    std::printf("File: %s, %d shells\n\n", AtomicTransitionManager::getEndfFileName("atom", kZ).c_str(),
                eadl.getNumberOfShells());

    bool ok = true;

    struct Ref { const char* shell; double eV; };
    const Ref refs[] = {{"K", 88005}, {"L1", 15861}, {"L2", 15200}, {"L3", 13035}, {"M1", 3851},
                        {"M2", 3554}, {"M3", 3066}, {"M4", 2586}, {"M5", 2484}};
    std::printf("%-5s %12s %12s %8s %8s %8s\n", "shell", "Booklet[keV]", "ENDF[keV]", "diff", "omega", "e-");
    for (const auto& r : refs) {
        const int i = findShell(eadl, r.shell);
        const double ref = r.eV * Units::eV;
        const double be = eadl.getBindingEnergy(i);
        const double relDiff = be / ref - 1.0;
        if (std::abs(relDiff) >= 0.005) {
            ok = false;
        }
        std::printf("%-5s %12.3f %12.3f %+7.2f%% %8.4f %8.1f\n", r.shell, ref / Units::keV, be / Units::keV,
                    100.0 * relDiff, eadl.getShell(i).totalRadiative, eadl.getShell(i).electrons);
    }

    double electrons = 0.0;
    for (int i = 0; i < eadl.getNumberOfShells(); ++i) {
        electrons += eadl.getShell(i).electrons;
    }
    if (std::abs(electrons - kZ) > 1.0e-6) {
        ok = false;
    }
    std::printf("sum of occupation numbers = %.1f (Z = %d)\n", electrons, kZ);

    const int k = findShell(eadl, "K");
    const double omegaRelDiff = eadl.getShell(k).totalRadiative / 0.963 - 1.0;
    if (std::abs(omegaRelDiff) >= 0.01) {
        ok = false;
    }
    std::printf("\nomega_K (sum of radiative K probabilities) = %.4f   Krause 1979: 0.963   (%+.2f%%)\n",
                eadl.getShell(k).totalRadiative, 100.0 * omegaRelDiff);

    double minSum = 1.0;
    double maxSum = 0.0;
    for (int i = 0; i < eadl.getNumberOfShells(); ++i) {
        if (eadl.hasTransitions(i)) {
            const double sum = eadl.getShell(i).totalRadiative + eadl.getShell(i).totalAuger;
            minSum = std::min(minSum, sum);
            maxSum = std::max(maxSum, sum);
        }
    }
    if (minSum <= 0.9 || maxSum > 1.0 + 1.0e-6) {
        ok = false;
    }
    std::printf("radiative + Auger branch sums (closed channels pruned): min = %.4f, max = %.6f\n", minSum, maxSum);

    std::printf("\nK lines [keV] vs Bearden 1967 (measured):\n");
    std::printf("%-12s %10s %10s %10s %10s %10s %8s\n", "line", "Bearden", "used", "diff", "file ETR", "diff", "prob");
    struct LineDef { const char* name; const char* origin; double bearden; };
    const LineDef lines[] = {{"Ka1 (K-L3)", "L3", 74.969}, {"Ka2 (K-L2)", "L2", 72.804},
                             {"Kb1 (K-M3)", "M3", 84.936}, {"Kb3 (K-M2)", "M2", 84.450}};
    for (const auto& l : lines) {
        const auto t = findLine(eadl, k, findShell(eadl, l.origin));
        const double ref = l.bearden * Units::keV;
        const double relDiff = t.energy / ref - 1.0;
        if (std::abs(relDiff) >= 0.001) {
            ok = false;
        }
        std::printf("%-12s %10.3f %10.3f %+9.2f%% %10.3f %+7.2f%% %8.4f\n", l.name, ref / Units::keV, t.energy / Units::keV,
                    100.0 * relDiff, t.fileEnergy / Units::keV, 100.0 * (t.fileEnergy / ref - 1.0), t.probability);
    }
    std::printf("used = B_K - B_o from EPICS2017 binding energies; file ETR = old EADL value.\n");

    long nAuger = 0;
    long nClosed = 0;
    double maxShift = 0.0;
    for (int i = 0; i < eadl.getNumberOfShells(); ++i) {
        for (const auto& t : eadl.getShell(i).auger) {
            ++nAuger;
            if (t.energy <= 0.0) {
                ++nClosed;
            }
            maxShift = std::max(maxShift, std::abs(t.energy - t.fileEnergy));
        }
    }
    if (nClosed > 0) {
        ok = false;
    }
    std::printf("Active Auger transitions: %ld (all E > 0), max |recomputed - file| = %.1f eV\n",
                nAuger, maxShift / Units::eV);
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 2  Subshell selection
// ==============================================================================
bool testShellSelection(int nSamples) {
    printHeader("TEST 2  Subshell selection: ENDF/B-VIII.0 vs the XCOM edge jumps");
    std::ifstream in(std::string(TOYSIM_DATA_DIR) + "/xcom_Z82.txt");
    if (!in) {
        throw std::runtime_error("cannot open xcom_Z82.txt");
    }

    bool ok = true;

    std::vector<std::pair<double, double>> rows; // (E [GeV], photoelectric mu/rho)
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        double e = 0.0, coh = 0.0, incoh = 0.0, pe = 0.0;
        if (s >> e >> coh >> incoh >> pe) {
            rows.emplace_back(e * Units::MeV, pe);
        }
    }

    const PhotoElectricShellData& pe = PhotoElectricShellData::get(kZ);
    const AtomicTransitionManager& atom = AtomicTransitionManager::get(kZ);
    std::printf("%d subshell cross sections; binding energies equal to the relaxation file: ", pe.getNumberOfShells());
    double worstB = 0.0;
    for (int i = 0; i < pe.getNumberOfShells(); ++i) {
        const int j = atom.findShellIndex(pe.getDesignator(i));
        if (j >= 0) {
            worstB = std::max(worstB, std::abs(pe.getBindingEnergy(i) / atom.getBindingEnergy(j) - 1.0));
        }
    }
    if (worstB >= 1.0e-6) {
        ok = false;
    }
    std::printf("max rel. diff %.1e\n\n", worstB);

    std::printf("%-5s %10s %10s | %10s %10s %8s\n", "shell", "XCOM edge", "1-lo/hi", "ENDF edge", "fraction", "diff");
    for (const char* name : {"K", "L1", "L2", "L3"}) {
        const int s = findPhotoShell(pe, name);
        const double edgeEpics = pe.getBindingEnergy(s);
        double best = 1.0e9, jump = 0.0, edgeXcom = 0.0;
        for (std::size_t i = 1; i < rows.size(); ++i) {
            if (rows[i].first == rows[i - 1].first && std::abs(rows[i].first - edgeEpics) < best) {
                best = std::abs(rows[i].first - edgeEpics);
                edgeXcom = rows[i].first;
                jump = 1.0 - rows[i - 1].second / rows[i].second;
            }
        }
        const double fraction = pe.computeShellFraction(s, edgeEpics * 1.0001);
        const double relDiff = fraction / jump - 1.0;
        if (std::abs(relDiff) >= 0.05) {
            ok = false;
        }
        std::printf("%-5s %10.3f %10.4f | %10.3f %10.4f %+7.2f%%\n", name, edgeXcom / Units::keV, jump,
                    edgeEpics / Units::keV, fraction, 100.0 * relDiff);
    }

    std::printf("\nSampled shell frequencies vs sigma_i / sum sigma (%d photons each):\n", nSamples);
    std::printf("%10s %6s %10s %10s %8s\n", "E [keV]", "shell", "computed", "sampled", "pull");
    for (double eKeV : {20.0, 100.0, 300.0, 1000.0}) {
        const double e = eKeV * Units::keV;
        std::vector<int> counts(static_cast<std::size_t>(pe.getNumberOfShells()), 0);
        for (int i = 0; i < nSamples; ++i) {
            const int s = pe.sampleShell(e);
            if (s >= 0) ++counts[static_cast<std::size_t>(s)];
        }
        for (const char* name : {"K", "L1", "L3", "M1"}) {
            const int s = findPhotoShell(pe, name);
            const double p = pe.computeShellFraction(s, e);
            if (p <= 0.0) continue;
            const double f = static_cast<double>(counts[static_cast<std::size_t>(s)]) / nSamples;
            const double pull = (f - p) / std::sqrt(p * (1.0 - p) / nSamples);
            if (std::abs(pull) >= kMaxPull) {
                ok = false;
            }
            std::printf("%10.0f %6s %10.4f %10.4f %+8.2f\n", eKeV, name, p, f, pull);
        }
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 3  Cascade sampling for K vacancies
// ==============================================================================
bool testCascade(int nVacancies) {
    printHeader("TEST 3  Relaxation cascade of a K vacancy (fluorescence + Auger)");
    const AtomicTransitionManager& atom = AtomicTransitionManager::get(kZ);
    const int k = findShell(atom, "K");
    AtomicDeexcitation deex(&gammaDef, &electronDef);
    deex.setAuger(true);

    struct Line {
        const char* name;
        double energy;
        double probability;
        long count = 0;
    };

    const auto makeLine = [&](const char* name, const char* originName) -> Line {
        const auto tr = findLine(atom, k, findShell(atom, originName));
        return {name, tr.energy, tr.probability, 0};
    };

    std::vector<Line> lines = {
        makeLine("Ka1 K-L3", "L3"),
        makeLine("Ka2 K-L2", "L2"),
        makeLine("Kb1 K-M3", "M3"),
        makeLine("Kb3 K-M2", "M2")
    };

    long kxrays = 0, lxrays = 0, augers = 0, photons = 0;
    std::vector<Track> secondaries;
    secondaries.reserve(16);

    for (int n = 0; n < nVacancies; ++n) {
        secondaries.clear();
        static_cast<void>(deex.generateParticles(kZ, k, Vector3D(), secondaries));
        for (const Track& t : secondaries) {
            if (t.getParticleDefinition()->getPdgId() == 11) {
                ++augers;
                continue;
            }
            ++photons;
            const double e = t.getKineticEnergy();
            if (e > 50.0 * Units::keV) {
                ++kxrays;
                for (auto& l : lines) {
                    if (std::abs(e - l.energy) < 1.0 * Units::eV) {
                        ++l.count;
                    }
                }
            } else if (e > 9.0 * Units::keV) {
                ++lxrays;
            }
        }
    }

    bool ok = true;
    const double omega = atom.getShell(k).totalRadiative;
    const double f = static_cast<double>(kxrays) / nVacancies;
    const double kPull = (f - omega) / std::sqrt(omega * (1.0 - omega) / nVacancies);
    if (std::abs(kPull) >= kMaxPull) {
        ok = false;
    }

    std::printf("%d K vacancies\n", nVacancies);
    std::printf("K X-rays per K vacancy: %.4f   omega_K = %.4f   pull %+.2f\n", f, omega, kPull);
    std::printf("%-10s %10s %10s %8s\n", "line", "data", "sampled", "pull");
    for (const auto& l : lines) {
        const double s = static_cast<double>(l.count) / nVacancies;
        const double pull = (s - l.probability) / std::sqrt(l.probability * (1.0 - l.probability) / nVacancies);
        if (std::abs(pull) >= kMaxPull) {
            ok = false;
        }
        std::printf("%-10s %10.4f %10.4f %+8.2f\n", l.name, l.probability, s, pull);
    }
    std::printf("Per K vacancy: %.3f photons (%.3f L X-rays 9-50 keV), %.3f Auger electrons\n",
                static_cast<double>(photons) / nVacancies, static_cast<double>(lxrays) / nVacancies,
                static_cast<double>(augers) / nVacancies);
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 4  Energy conservation
// ==============================================================================
bool testEnergyConservation(int nEvents) {
    printHeader("TEST 4  Energy conservation: E = T_e + X-rays + Auger + local (event by event)");
    auto lead = makeLead();
    std::printf("%-26s %8s %12s %14s %10s %10s\n", "setup", "E [keV]", "max|dE|/E", "min local", "sec/event", "local/E");

    AtomicDeexcitation deex(&gammaDef, &electronDef);
    PhotoElectricEffect phot(&electronDef);
    phot.setVerboseLevel(0);
    phot.setAtomicDeexcitation(&deex);

    std::vector<Track> secondaries;
    secondaries.reserve(16);

    double overallWorst = 0.0;
    double overallMinLocal = 1.0e9;

    for (bool augerOn : {true, false}) {
        deex.setAuger(augerOn);
        for (double eKeV : {20.0, 89.0, 100.0, 1000.0}) {
            const double e = eKeV * Units::keV;
            double worst = 0.0, minLocal = 1.0e9, sumLocal = 0.0;
            long nSec = 0;
            for (int n = 0; n < nEvents; ++n) {
                Track photon(&gammaDef, Vector3D(), Vector3D(0.0, 0.0, 1.0), e);
                Step step;
                secondaries.clear();
                phot.postStepDoIt(photon, *lead, step, secondaries);

                double sum = step.totalEnergyDeposit;
                for (const Track& t : secondaries) {
                    sum += t.getKineticEnergy();
                }
                worst = std::max(worst, std::abs(sum - e) / e);
                minLocal = std::min(minLocal, step.totalEnergyDeposit);
                sumLocal += step.totalEnergyDeposit;
                nSec += static_cast<long>(secondaries.size());
            }
            overallWorst = std::max(overallWorst, worst);
            overallMinLocal = std::min(overallMinLocal, minLocal);

            char label[64];
            std::snprintf(label, sizeof(label), "ENDF/B-VIII.0, Auger %s", augerOn ? "on" : "off");
            std::printf("%-26s %8.0f %12.2e %10.3f keV %10.2f %10.4f\n", label, eKeV, worst, minLocal / Units::keV,
                        static_cast<double>(nSec) / nEvents, sumLocal / (nEvents * e));
        }
    }

    const bool ok = (overallWorst < 1.0e-12) && (overallMinLocal >= 0.0);
    if (overallMinLocal >= 0.0) {
        std::printf("min local >= 0: no transition emits more than the binding energy available.\n");
    } else {
        std::printf("FAIL: negative local energy deposit encountered (min local = %.6f keV).\n",
                    overallMinLocal / Units::keV);
    }
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 5  K X-ray escape from a semi-infinite lead block
// ==============================================================================
bool testEscape(int nPhotons) {
    printHeader("TEST 5  K X-ray escape from a thick lead block vs the analytic formula");
    auto lead = makeLead();
    const AtomicTransitionManager& atom = AtomicTransitionManager::get(kZ);
    const PhotoElectricShellData& pe = PhotoElectricShellData::get(kZ);
    const int k = findShell(atom, "K");
    const int kPhoto = findPhotoShell(pe, "K");
    PhotoElectricEffect probe(&electronDef);

    bool ok = true;
    std::printf("%8s %8s %10s %10s %10s %8s %10s\n", "E0[keV]", "f_K", "MC", "formula", "diff", "pull", "L esc");
    for (double eKeV : {100.0, 150.0, 300.0}) {
        const double e0 = eKeV * Units::keV;
        const double mu0 = probe.computeMassAttenuation(e0, *lead) * lead->getDensity();
        const double fK = pe.computeShellFraction(kPhoto, e0);
        double expected = 0.0;
        for (const auto& t : atom.getShell(k).radiative) {
            const double mui = probe.computeMassAttenuation(t.energy, *lead) * lead->getDensity();
            expected += 0.5 * fK * t.probability * (1.0 - (mui / mu0) * std::log1p(mu0 / mui));
        }

        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        rm.setWorldVolume(std::make_unique<Box>("World", 31.0, 31.0, 11.0, vacuum));
        rm.addVolume(std::make_unique<Box>("Lead", 30.0, 30.0, 10.0, lead));
        rm.addPhysicsProcess(std::make_unique<Ionization>());

        AtomicDeexcitation deex(&gammaDef, &electronDef);
        auto phot = std::make_unique<PhotoElectricEffect>(&electronDef);
        phot->setAtomicDeexcitation(&deex);
        rm.addPhysicsProcess(std::move(phot));

        rm.setEnergyCut(1.0 * Units::keV);
        rm.setVerboseLevel(0);

        long kEscaped = 0, lEscaped = 0;
        rm.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (!step.leftWorld || track.getParticleDefinition()->getPdgId() != 22 || track.getDirection().z >= 0.0) {
                return;
            }
            const double e = track.getKineticEnergy();
            if (e > 60.0 * Units::keV && e < 88.0 * Units::keV) ++kEscaped;
            else if (e > 9.0 * Units::keV && e < 16.0 * Units::keV) ++lEscaped;
        });

        for (int i = 0; i < nPhotons; ++i) {
            Track photon(&gammaDef, Vector3D(0.0, 0.0, -10.5), Vector3D(0.0, 0.0, 1.0), e0);
            rm.simulateEvent(photon, 1.0 * Units::m, 100000);
        }

        const double mc = static_cast<double>(kEscaped) / nPhotons;
        const double sigma = std::sqrt(expected * (1.0 - expected) / nPhotons);
        const double pull = (mc - expected) / sigma;
        if (std::abs(pull) >= kMaxPull) {
            ok = false;
        }
        std::printf("%8.0f %8.4f %10.5f %10.5f %+9.2f%% %+8.2f %10.5f\n", eKeV, fK, mc, expected,
                    100.0 * (mc / expected - 1.0), pull, static_cast<double>(lEscaped) / nPhotons);
    }
    std::printf("Note: photoelectric-only (no Compton or Rayleigh scattering, matching the formula's mu_0 and mu_i).\n");
    std::printf("MC / formula: K X-rays (60-88 keV) leaving through the front face per incident photon.\n");
    std::printf("L esc: L X-rays (9-16 keV) per incident photon (no formula: they also come from the cascade).\n");
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

} // namespace

int main(int argc, char** argv) {
    bool quick = false;
    bool escapeOnly = false;
    unsigned long long seed = 12345;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "quick") quick = true;
        else if (arg == "escape") escapeOnly = true;
        else seed = std::stoull(arg);
    }
    RandomEngine::setSeed(seed);

    std::printf("[RandomEngine] Seed: %llu%s\n", static_cast<unsigned long long>(RandomEngine::getSeed()),
                quick ? "   (quick mode)" : "");

    if (escapeOnly) {
        const bool ok = testEscape(quick ? 5000 : 100000);
        std::printf("\n==============================================================================\n");
        std::printf("  OVERALL: %s\n", ok ? "PASS" : "FAIL");
        std::printf("==============================================================================\n");
        return ok ? 0 : 1;
    }

    bool allPassed = true;
    allPassed &= testShellData();
    allPassed &= testShellSelection(quick ? 50000 : 1000000);
    allPassed &= testCascade(quick ? 50000 : 1000000);
    allPassed &= testEnergyConservation(quick ? 2000 : 20000);
    allPassed &= testEscape(quick ? 5000 : 100000);

    std::printf("\n==============================================================================\n");
    std::printf("  OVERALL: %s\n", allPassed ? "PASS" : "FAIL");
    std::printf("==============================================================================\n");
    return allPassed ? 0 : 1;
}