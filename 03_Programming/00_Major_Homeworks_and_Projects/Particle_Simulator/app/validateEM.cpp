// ==============================================================================
// [STAGE 3] Validation of the photon and electron physics, and of the shower.
//
//   TEST 1  Photon cross sections (mu/rho) vs NIST XCOM
//   TEST 2  Narrow-beam attenuation: Monte Carlo vs exp(-mu x)
//   TEST 3  Sampled distributions vs the differential cross sections
//   TEST 4  Kinematics of one interaction: energy and momentum conservation
//   TEST 5  Electrons: radiative / total stopping power and CSDA range vs NIST ESTAR
//   TEST 6  Energy conservation in a full shower
//   TEST 7  Longitudinal shower profile vs the PDG parametrisation
//
// Run:  ./ValidationEMApp            (all tests, ~1 min)
//       ./ValidationEMApp quick      (fewer events, a few seconds)
// ==============================================================================
#include "Annihilation.h"
#include "Box.h"
#include "Bremsstrahlung.h"
#include "ComptonScattering.h"
#include "GammaConversion.h"
#include "Ionization.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "Particle.h"
#include "PhotoElectricEffect.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

using PhysicalConstants::m_e;

namespace {

// ------------------------------------------------------------------------------
// Particles and the material
// ------------------------------------------------------------------------------
const Particle gammaDef(22, "gamma", 0.0, 0.0, 1.0, -1.0);
const Particle electronDef(11, "e-", m_e, -1.0, 0.5, -1.0);
const Particle positronDef(-11, "e+", m_e, 1.0, 0.5, -1.0);

std::shared_ptr<Material> makeLead() {
    DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

// cm^2 per atom -> cm^2/g
double perGram(double sigmaPerAtom, const Material& material) {
    return sigmaPerAtom * PhysicalConstants::N_A / material.getA();
}

void header(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

// ------------------------------------------------------------------------------
// Data tables
// ------------------------------------------------------------------------------
struct XcomRow { double energy, coherent, incoherent, photoelectric, pair, total; };  // energy in GeV, rest cm2/g
struct EstarRow { double energy, collision, radiative, total, range; };             // MeV cm2/g, g/cm2

std::vector<std::vector<double>> readColumns(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream s(line);
        std::vector<double> row;
        double v;
        while (s >> v) row.push_back(v);
        if (!row.empty()) rows.push_back(row);
    }
    return rows;
}

std::vector<XcomRow> readXcom() {
    std::vector<XcomRow> rows;
    for (const auto& r : readColumns(std::string(TOYSIM_DATA_DIR) + "/xcom_Z82.txt")) {
        rows.push_back({r[0] * Units::MeV, r[1], r[2], r[3], r[4], r[5]});
    }
    return rows;
}

std::vector<EstarRow> readEstar() {
    std::vector<EstarRow> rows;
    for (const auto& r : readColumns(std::string(TOYSIM_DATA_DIR) + "/estar_Pb.txt")) {
        rows.push_back({r[0] * Units::MeV, r[1], r[2], r[3], r[4]});
    }
    return rows;
}

// ------------------------------------------------------------------------------
// The full EM physics list (like G4EmStandardPhysics, for e-, e+, gamma)
// ------------------------------------------------------------------------------
struct EmSetup {
    RunManager runManager;
    std::shared_ptr<Material> lead;
};

// A lead box of half-sizes (hx, hy, hz) at the origin, inside a vacuum world
// 1 cm larger on every side.
std::unique_ptr<EmSetup> makeSetup(double hx, double hy, double hz, double energyCut = 100.0 * Units::keV) {
    auto setup = std::make_unique<EmSetup>();
    setup->lead = makeLead();
    auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
    RunManager& rm = setup->runManager;
    rm.setWorldVolume(std::make_unique<Box>("World", hx + 1.0, hy + 1.0, hz + 1.0, vacuum));
    rm.addVolume(std::make_unique<Box>("Lead", hx, hy, hz, setup->lead));

    // msc first (pre-step momentum), then the continuous losses, then the discrete ones
    rm.addPhysicsProcess(std::make_unique<MultipleScattering>());
    rm.addPhysicsProcess(std::make_unique<Ionization>());
    rm.addPhysicsProcess(std::make_unique<Bremsstrahlung>(&gammaDef));
    rm.addPhysicsProcess(std::make_unique<Annihilation>(&gammaDef));
    rm.addPhysicsProcess(std::make_unique<PhotoElectricEffect>(&electronDef));
    rm.addPhysicsProcess(std::make_unique<ComptonScattering>(&electronDef));
    rm.addPhysicsProcess(std::make_unique<GammaConversion>(&electronDef, &positronDef));
    rm.setEnergyCut(energyCut);
    rm.setVerboseLevel(0);
    return setup;
}

// chi^2 of a histogram against expected counts (bins with expected < 5 skipped)
void printChi2(const char* label, const std::vector<double>& observed, const std::vector<double>& expected) {
    double chi2 = 0.0;
    int ndf = 0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        if (expected[i] < 5.0) continue;
        chi2 += (observed[i] - expected[i]) * (observed[i] - expected[i]) / expected[i];
        ++ndf;
    }
    std::printf("  %-36s chi2/ndf = %7.1f / %3d = %5.2f\n", label, chi2, ndf, chi2 / std::max(ndf, 1));
}

// ==============================================================================
// TEST 1  mu/rho vs XCOM
// ==============================================================================
void testCrossSections(const std::vector<XcomRow>& xcom) {
    header("TEST 1  Photon mass attenuation coefficients (cm^2/g) vs NIST XCOM, lead");
    std::printf("photoelectric and pair: the simulation USES the XCOM table (diff = interpolation only).\n");
    std::printf("Compton: Klein-Nishina for free electrons.  BH: pair cross section from the integrated\n");
    std::printf("Bethe-Heitler formula (the one that samples eps), for comparison only.\n");
    std::printf("Total = without coherent (Rayleigh) scattering, which we do not simulate.\n\n");
    std::printf("%9s | %9s %7s | %9s %9s %7s | %9s %7s | %9s %7s | %9s %9s %7s\n",
                "E [MeV]", "phot", "diff", "Compton", "XCOM", "diff", "pair", "diff", "BH", "diff", "total", "XCOM", "diff");

    auto lead = makeLead();
    PhotoElectricEffect phot(&electronDef);
    phot.setVerboseLevel(0);
    const double energies[] = {0.01, 0.03, 0.1, 0.3, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0, 1000.0, 10000.0, 100000.0};
    for (double eMeV : energies) {
        const double energy = eMeV * Units::MeV;
        const auto row = std::find_if(xcom.begin(), xcom.end(), [&](const XcomRow& r) { return std::abs(r.energy / energy - 1.0) < 1e-6; });
        if (row == xcom.end()) continue;
        const double pe = phot.computeMassAttenuation(energy, *lead);
        const double co = perGram(lead->getZ() * ComptonScattering::computeCrossSectionPerElectron(energy), *lead);
        const double pa = perGram(GammaConversion::getCrossSectionPerAtom(energy, *lead), *lead);
        const double bh = perGram(GammaConversion::computeCrossSectionPerAtom(energy, *lead), *lead);
        const double tot = pe + co + pa;
        auto diff = [](double a, double b) { return b > 0.0 ? 100.0 * (a / b - 1.0) : 0.0; };
        std::printf("%9.3g | %9.3e %+6.1f%% | %9.3e %9.3e %+6.1f%% | %9.3e %+6.1f%% | %9.3e %+6.1f%% | %9.3e %9.3e %+6.1f%%\n",
                    eMeV, pe, diff(pe, row->photoelectric), co, row->incoherent, diff(co, row->incoherent),
                    pa, diff(pa, row->pair), bh, diff(bh, row->pair), tot, row->total, diff(tot, row->total));
    }
    std::printf("\nHigh-energy limit of pair production: lambda -> 9/7 X0 = 1.286 X0  (PDG X0 = 6.37 g/cm^2)\n");
    for (double eGeV : {1.0, 10.0, 100.0}) {
        const double mu = perGram(GammaConversion::getCrossSectionPerAtom(eGeV * Units::GeV, *lead), *lead);
        const double muBH = perGram(GammaConversion::computeCrossSectionPerAtom(eGeV * Units::GeV, *lead), *lead);
        std::printf("  E = %4.0f GeV : lambda_pair = %.3f X0 (XCOM), %.3f X0 (Bethe-Heitler)\n", eGeV, 1.0 / mu / 6.37, 1.0 / muBH / 6.37);
    }
}

// ==============================================================================
// TEST 2  Narrow-beam attenuation
// ==============================================================================
// A pencil beam on a slab of thickness x. A photon that crosses without ANY
// interaction leaves with its energy and direction unchanged. Their fraction
// must be exp(-mu x), mu = the sum of our three processes. This checks the whole
// chain: the three GPILs competing, the "number of mean free paths" bookkeeping
// across steps, and the geometry.
void testAttenuation(int nPhotons) {
    header("TEST 2  Narrow-beam attenuation through a lead slab: MC vs exp(-mu x)");
    std::printf("%9s %8s %8s | %10s %10s %7s %8s\n", "E [MeV]", "x [cm]", "mu x", "MC", "exp(-mu x)", "diff", "pull");

    for (double eMeV : {0.1, 0.5, 2.0, 10.0, 100.0}) {
        const double energy = eMeV * Units::MeV;
        auto lead = makeLead();
        PhotoElectricEffect phot(&electronDef);
        phot.setVerboseLevel(0);
        const double muRho = phot.computeMassAttenuation(energy, *lead)
                           + perGram(lead->getZ() * ComptonScattering::computeCrossSectionPerElectron(energy), *lead)
                           + perGram(GammaConversion::getCrossSectionPerAtom(energy, *lead), *lead);
        const double mu = muRho * lead->getDensity();
        const double thickness = 1.5 / mu; // 1.5 mean free paths: ~22 % survive

        auto setup = makeSetup(50.0, 50.0, 0.5 * thickness);
        int uncollided = 0;
        setup->runManager.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (step.leftWorld && track.getParticleDefinition()->getPdgId() == 22
                && track.getKineticEnergy() == energy && track.getDirection().z == 1.0) {
                ++uncollided;
            }
        });
        for (int i = 0; i < nPhotons; ++i) {
            Track photon(&gammaDef, Vector3D(0, 0, -0.5 * thickness - 0.5), Vector3D(0, 0, 1), energy);
            setup->runManager.simulateEvent(photon, 1.0 * Units::m, 100000);
        }
        const double expected = std::exp(-mu * thickness);
        const double measured = static_cast<double>(uncollided) / nPhotons;
        const double sigma = std::sqrt(expected * (1.0 - expected) / nPhotons);
        std::printf("%9.3g %8.4f %8.3f | %10.4f %10.4f %+6.2f%% %+8.2f\n", eMeV, thickness, mu * thickness,
                    measured, expected, 100.0 * (measured / expected - 1.0), (measured - expected) / sigma);
    }
    std::printf("pull = (MC - expected) / binomial sigma: should be within ~ +-2.\n");
}

// ==============================================================================
// TEST 3  Sampled distributions vs the formulas they come from
// ==============================================================================
// Each sampler is called N times directly; the histogram is compared bin by bin
// with N * (integral of the normalised differential cross section over the bin).
template <typename Pdf>
std::vector<double> expectedCounts(Pdf pdf, double lo, double hi, int nBins, int nSamples, bool logBins) {
    std::vector<double> integral(nBins, 0.0);
    double total = 0.0;
    for (int b = 0; b < nBins; ++b) {
        const double a = logBins ? lo * std::pow(hi / lo, static_cast<double>(b) / nBins) : lo + (hi - lo) * b / nBins;
        const double c = logBins ? lo * std::pow(hi / lo, static_cast<double>(b + 1) / nBins) : lo + (hi - lo) * (b + 1) / nBins;
        const int n = 40;
        double s = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double x = a + (c - a) * i / n;
            s += ((i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0)) * pdf(x);
        }
        integral[b] = s * (c - a) / (3.0 * n);
        total += integral[b];
    }
    for (double& v : integral) v *= nSamples / total;
    return integral;
}

void fill(std::vector<double>& histogram, double x, double lo, double hi, bool logBins) {
    const int nBins = static_cast<int>(histogram.size());
    const double f = logBins ? std::log(x / lo) / std::log(hi / lo) : (x - lo) / (hi - lo);
    const int b = static_cast<int>(f * nBins);
    if (b >= 0 && b < nBins) histogram[b] += 1.0;
}

void testSampling(int nSamples) {
    header("TEST 3  Sampled distributions vs the differential cross sections");
    auto lead = makeLead();
    const int nBins = 50;

    std::printf("Compton, eps = E'/E (Klein-Nishina):\n");
    for (double eMeV : {0.1, 1.0, 10.0}) {
        const double energy = eMeV * Units::MeV;
        const double lo = 1.0 / (1.0 + 2.0 * energy / m_e), hi = 1.0;
        std::vector<double> h(nBins, 0.0);
        for (int i = 0; i < nSamples; ++i) fill(h, ComptonScattering::sampleEnergyFraction(energy), lo, hi, false);
        const auto e = expectedCounts([&](double x) { return ComptonScattering::computeDifferentialCrossSection(energy, x); },
                                      lo, hi, nBins, nSamples, false);
        char label[64];
        std::snprintf(label, sizeof(label), "E = %g MeV", eMeV);
        printChi2(label, h, e);
    }

    std::printf("Pair production, eps = E(e)/E_gamma folded to [eps0, 1/2] (Bethe-Heitler):\n");
    for (double eMeV : {5.0, 100.0, 10000.0}) {
        const double energy = eMeV * Units::MeV;
        const double lo = m_e / energy, hi = 0.5;
        std::vector<double> h(nBins, 0.0);
        for (int i = 0; i < nSamples; ++i) fill(h, GammaConversion::sampleEnergyFraction(energy, *lead), lo, hi, false);
        const auto e = expectedCounts([&](double x) { return GammaConversion::computeDifferentialCrossSection(energy, x, *lead); },
                                      lo, hi, nBins, nSamples, false);
        char label[64];
        std::snprintf(label, sizeof(label), "E = %g MeV", eMeV);
        printChi2(label, h, e);
    }

    std::printf("Bremsstrahlung, photon energy k in [k_c = 100 keV, T] (log bins):\n");
    for (double tMeV : {10.0, 1000.0}) {
        const double t = tMeV * Units::MeV, kc = 100.0 * Units::keV;
        std::vector<double> h(nBins, 0.0);
        for (int i = 0; i < nSamples; ++i) fill(h, Bremsstrahlung::samplePhotonEnergy(t, kc, *lead), kc, t, true);
        const auto e = expectedCounts([&](double k) { return Bremsstrahlung::computeDifferentialCrossSection(t, k, *lead); },
                                      kc, t, nBins, nSamples, true);
        char label[64];
        std::snprintf(label, sizeof(label), "T = %g MeV", tMeV);
        printChi2(label, h, e);
    }
    std::printf("chi2/ndf ~ 1 means the sampler reproduces its formula.\n");
}

// ==============================================================================
// TEST 4  One interaction at a time: conservation laws
// ==============================================================================
// Call the postStepDoIt of each process many times and check, for every single
// interaction, that energy (and where the target recoil is negligible, momentum)
// is conserved. Errors here are kinematics bugs: a wrong cos(theta) formula,
// a sign, a missing mass.
void testKinematics(int nSamples) {
    header("TEST 4  Kinematics of single interactions (max violation over all samples)");
    auto lead = makeLead();
    std::printf("%-32s %12s %12s\n", "", "max |dE|/E", "max |dp|/p");

    // targetMass: rest energy of the target that is used up or set free. Compton and
    // photoelectric eject an electron that was already there (its m_e existed before);
    // annihilation destroys an atomic electron (m_e set free as photon energy).
    auto run = [&](const char* label, PhysicsProcess& process, const Particle& particle, double energy,
                   bool checkMomentum, double targetMass) {
        process.setVerboseLevel(0);
        double maxDE = 0.0, maxDP = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            Track track(&particle, Vector3D(), Vector3D(0.3, -0.2, 1.0), energy);
            const double eBefore = track.getTotalEnergy() + targetMass;
            const Vector3D pBefore = track.getMomentum();
            Step step;
            std::vector<Track> secondaries;
            process.postStepDoIt(track, *lead, step, secondaries);
            double eAfter = step.totalEnergyDeposit + (track.isAlive() ? track.getTotalEnergy() : 0.0);
            Vector3D pAfter = track.isAlive() ? track.getMomentum() : Vector3D();
            for (const auto& s : secondaries) {
                eAfter += s.getTotalEnergy();
                pAfter = pAfter + s.getMomentum();
            }
            maxDE = std::max(maxDE, std::abs(eAfter - eBefore) / eBefore);
            maxDP = std::max(maxDP, (pAfter - pBefore).mag() / pBefore.mag());
        }
        if (checkMomentum) std::printf("%-32s %12.2e %12.2e\n", label, maxDE, maxDP);
        else std::printf("%-32s %12.2e %12s\n", label, maxDE, "(nucleus)");
    };

    ComptonScattering compt(&electronDef);
    PhotoElectricEffect phot(&electronDef);
    GammaConversion conv(&electronDef, &positronDef);
    Bremsstrahlung brem(&gammaDef);
    Annihilation annihil(&gammaDef);
    run("Compton 1 MeV", compt, gammaDef, 1.0 * Units::MeV, true, m_e);
    run("Compton 100 keV", compt, gammaDef, 0.1 * Units::MeV, true, m_e);
    run("photoelectric 200 keV", phot, gammaDef, 0.2 * Units::MeV, false, m_e);
    run("pair production 10 MeV", conv, gammaDef, 10.0 * Units::MeV, false, 0.0);
    run("bremsstrahlung 1 GeV", brem, electronDef, 1.0 * Units::GeV, false, 0.0);
    run("annihilation in flight 10 MeV", annihil, positronDef, 10.0 * Units::MeV, true, m_e);
    std::printf("Compton and annihilation are two-body: momentum exact. Photoelectric, pair\n");
    std::printf("and brems give momentum to the nucleus (not simulated): energy only.\n");
}

// ==============================================================================
// TEST 5  Electrons vs ESTAR
// ==============================================================================
void testElectrons(const std::vector<EstarRow>& estar) {
    header("TEST 5  Electrons in lead vs NIST ESTAR (MeV cm^2/g, g/cm^2)");
    auto lead = makeLead();
    Ionization ion;
    ion.setVerboseLevel(0);
    const double rho = lead->getDensity();
    auto collision = [&](double t) { return ion.computeDEDX(electronDef, t, *lead) / rho / Units::MeV; };
    auto radiative = [&](double t) { return Bremsstrahlung::computeRadiativeDEDX(t, *lead) / rho / Units::MeV; };

    // CSDA range with the TOTAL loss: R(T) = R_col(1 keV) + integral_{1 keV}^{T} dT / (S_col + S_rad)
    // (radiative loss is negligible below 1 keV). Trapezoid rule on a fine ln(T) grid.
    std::vector<double> gridT, gridR;
    {
        const double t0 = 1.0 * Units::keV, t1 = 1.0 * Units::GeV;
        const int n = 4000;
        const double h = std::log(t1 / t0) / n;
        double range = ion.getRange(electronDef, t0, *lead) * rho;
        double prev = t0 / ((collision(t0) + radiative(t0)) * Units::MeV);
        gridT.push_back(t0);
        gridR.push_back(range);
        for (int i = 1; i <= n; ++i) {
            const double t = t0 * std::exp(i * h);
            const double cur = t / ((collision(t) + radiative(t)) * Units::MeV);
            range += 0.5 * (prev + cur) * h;
            prev = cur;
            gridT.push_back(t);
            gridR.push_back(range);
        }
    }
    auto csdaRange = [&](double t) {
        const auto it = std::lower_bound(gridT.begin(), gridT.end(), t);
        const std::size_t j = std::max<std::size_t>(1, it - gridT.begin());
        const double f = std::log(t / gridT[j - 1]) / std::log(gridT[j] / gridT[j - 1]);
        return gridR[j - 1] + f * (gridR[j] - gridR[j - 1]);
    };

    std::printf("%9s | %9s %9s %7s | %9s %9s %7s | %9s %9s %7s\n", "T [MeV]", "S_rad", "ESTAR", "diff",
                "S_total", "ESTAR", "diff", "R_CSDA", "ESTAR", "diff");
    for (double tMeV : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0, 1000.0}) {
        const auto row = std::find_if(estar.begin(), estar.end(), [&](const EstarRow& r) { return std::abs(r.energy / (tMeV * Units::MeV) - 1.0) < 1e-6; });
        if (row == estar.end()) continue;
        const double t = row->energy;
        const double sr = radiative(t), st = collision(t) + sr, r = csdaRange(t);
        std::printf("%9.3g | %9.3e %9.3e %+6.1f%% | %9.3e %9.3e %+6.1f%% | %9.3e %9.3e %+6.1f%%\n", tMeV,
                    sr, row->radiative, 100.0 * (sr / row->radiative - 1.0),
                    st, row->total, 100.0 * (st / row->total - 1.0),
                    r, row->range, 100.0 * (r / row->range - 1.0));
    }

    // Critical energy, two definitions (PDG section 34.4.4):
    //   S_col(Ec) = S_rad(Ec)                              (Berger & Seltzer, ESTAR)
    //   S_col(Ec) = Ec / X0   (Rossi: the loss per X0 equals the energy)   PDG: 7.43 MeV for lead
    auto bisect = [](std::function<double(double)> f, double lo, double hi) {
        for (int i = 0; i < 100; ++i) {
            const double mid = std::sqrt(lo * hi);
            (f(lo) * f(mid) <= 0.0 ? hi : lo) = mid;
        }
        return std::sqrt(lo * hi);
    };
    const double x0 = 6.37; // g/cm^2
    const double ecEqual = bisect([&](double t) { return collision(t) - radiative(t); }, 1.0 * Units::MeV, 100.0 * Units::MeV);
    const double ecRossi = bisect([&](double t) { return collision(t) - t / Units::MeV / x0; }, 1.0 * Units::MeV, 100.0 * Units::MeV);
    double ecEstar = 0.0;
    for (std::size_t i = 1; i < estar.size(); ++i) {
        const double d0 = estar[i - 1].collision - estar[i - 1].radiative, d1 = estar[i].collision - estar[i].radiative;
        if (d0 > 0.0 && d1 <= 0.0) ecEstar = estar[i - 1].energy + (estar[i].energy - estar[i - 1].energy) * d0 / (d0 - d1);
    }
    std::printf("\nCritical energy, S_col = S_rad      : ours %.2f MeV, ESTAR %.2f MeV\n", ecEqual / Units::MeV, ecEstar / Units::MeV);
    std::printf("Critical energy, Rossi (S_col = E/X0): ours %.2f MeV, PDG 7.43 MeV\n", ecRossi / Units::MeV);
    std::printf("High-energy brems, S_rad / (E/X0) with PDG X0 = 6.37 g/cm^2:\n ");
    for (double tGeV : {0.1, 1.0, 10.0, 100.0}) {
        const double t = tGeV * Units::GeV;
        std::printf("  %g GeV: %.3f", tGeV, radiative(t) * Units::MeV / ((t + m_e) / x0));
    }
    std::printf("\n (-> 1 + 1/(18 ln(183 Z^-1/3)) = 1.015 at full screening: the term X0 leaves out)\n");
}

// ==============================================================================
// Shower bookkeeping shared by tests 6 and 7
// ==============================================================================
struct ShowerTally {
    double deposited = 0.0;      // sum of step deposits
    double escaped = 0.0;        // kinetic energy of the particles that left the world
    int conversions = 0;         // pair productions: each puts 2 m_e into mass
    int annihilations = 0;       // each releases 2 m_e (e+ and an atomic e-)
    std::vector<double> profile; // deposit in bins of depth t = z / X0
    double binWidth = 0.25;      // in X0
    double zFront = 0.0;         // z of the lead's front face
};

void attachTally(RunManager& rm, ShowerTally& tally, double x0) {
    rm.setUserSteppingAction([&tally, x0](const Track& track, const Step& step) {
        tally.deposited += step.totalEnergyDeposit;
        if (step.leftWorld) tally.escaped += track.getKineticEnergy();
        if (step.processDefinedStep == "conv") ++tally.conversions;
        if (step.processDefinedStep == "annihil") ++tally.annihilations;
        if (step.totalEnergyDeposit > 0.0 && !tally.profile.empty()) {
            // deposit at the step's midpoint
            const double z = 0.5 * (step.preStepPosition.z + step.postStepPosition.z);
            const int b = static_cast<int>((z - tally.zFront) / x0 / tally.binWidth);
            if (b >= 0 && b < static_cast<int>(tally.profile.size())) tally.profile[b] += step.totalEnergyDeposit;
        }
    });
}

// ==============================================================================
// TEST 6  Energy conservation in a shower
// ==============================================================================
// For an incident particle of kinetic energy E0, everything must be accounted for:
//   deposited + escaped = E0 + 2 m_e (N_annihilations - N_conversions)
// exactly (rounding only). Each pair production turns 2 m_e of photon energy into
// mass; each annihilation turns 2 m_e of mass back into photons.
void testEnergyConservation(int nEvents) {
    header("TEST 6  Energy conservation in a full shower (lead block 20 x 20 x 20 cm)");
    std::printf("%-13s %6s | %10s %10s %8s %8s | %12s %8s %9s\n", "primary", "events", "deposited", "escaped",
                "N_conv", "N_annih", "max |imbal.|", "tracks", "max stack");
    struct Case { const Particle* p; double energy; const char* name; };
    for (const Case& c : {Case{&electronDef, 1.0 * Units::GeV, "e- 1 GeV"}, Case{&gammaDef, 1.0 * Units::GeV, "gamma 1 GeV"},
                          Case{&positronDef, 10.0 * Units::MeV, "e+ 10 MeV"}, Case{&gammaDef, 0.5 * Units::MeV, "gamma 0.5 MeV"}}) {
        auto setup = makeSetup(10.0, 10.0, 10.0);
        double sumDep = 0.0, sumEsc = 0.0, maxImbalance = 0.0;
        long sumConv = 0, sumAnnih = 0, sumTracks = 0;
        std::size_t maxStack = 0;
        ShowerTally tally;
        attachTally(setup->runManager, tally, setup->lead->getRadLength());
        for (int i = 0; i < nEvents; ++i) {
            tally = ShowerTally();
            Track primary(c.p, Vector3D(0, 0, -10.5), Vector3D(0, 0, 1), c.energy);
            setup->runManager.simulateEvent(primary, 1.0 * Units::m, 1000000);
            const double balance = tally.deposited + tally.escaped
                                 - (c.energy + 2.0 * m_e * (tally.annihilations - tally.conversions));
            maxImbalance = std::max(maxImbalance, std::abs(balance) / c.energy);
            sumDep += tally.deposited;
            sumEsc += tally.escaped;
            sumConv += tally.conversions;
            sumAnnih += tally.annihilations;
            sumTracks += setup->runManager.getNumberOfTracks();
            maxStack = std::max(maxStack, setup->runManager.getMaxStackSize());
        }
        std::printf("%-13s %6d | %8.3f %% %8.3f %% %8.1f %8.1f | %12.1e %8.0f %9zu\n", c.name, nEvents,
                    100.0 * sumDep / (nEvents * c.energy), 100.0 * sumEsc / (nEvents * c.energy),
                    static_cast<double>(sumConv) / nEvents, static_cast<double>(sumAnnih) / nEvents,
                    maxImbalance, static_cast<double>(sumTracks) / nEvents, maxStack);
    }
    std::printf("deposited/escaped: mean fraction of E0. imbalance: rounding only (~1e-15).\n");
    std::printf("tracks: mean number per event; max stack: largest track stack seen.\n");
}

// ==============================================================================
// TEST 7  Longitudinal profile
// ==============================================================================
// PDG eq. 34.36: dE/dt = E0 b (b t)^(a-1) e^(-b t) / Gamma(a), t = depth in X0,
// with  t_max = (a - 1)/b = ln(E0/Ec) + C,  C = -0.5 (e-), +0.5 (gamma), Ec = 7.43 MeV,
// and b ~ 0.5. We take the moments of OUR profile: <t> = a/b, var(t) = a/b^2 -> a, b.
void testShowerProfile(int nEvents1GeV) {
    header("TEST 7  Longitudinal shower profile in lead vs the PDG parametrisation");
    std::printf("t_max from the moments ((a-1)/b) and from the peak of the profile (3-bin smoothed).\n");
    std::printf("%-13s %6s | %6s %6s | %7s %7s %7s | %9s %9s %9s %7s\n", "primary", "events", "a", "b",
                "moments", "peak", "PDG", "contained", "tracks", "max stack", "s/evt");
    const double ec = 7.43 * Units::MeV;
    struct Case { const Particle* p; double energy; const char* name; double c; int events; };
    const std::vector<Case> cases = {
        {&electronDef, 1.0 * Units::GeV, "e- 1 GeV", -0.5, nEvents1GeV},
        {&electronDef, 10.0 * Units::GeV, "e- 10 GeV", -0.5, std::max(5, nEvents1GeV / 5)},
        {&gammaDef, 10.0 * Units::GeV, "gamma 10 GeV", 0.5, std::max(5, nEvents1GeV / 5)},
        {&electronDef, 100.0 * Units::GeV, "e- 100 GeV", -0.5, std::max(2, nEvents1GeV / 50)}};

    std::vector<double> shape10;
    for (const Case& c : cases) {
        // 40 X0 deep (22.4 cm), 1 m wide: practically full containment
        const double halfDepth = 20.0 * 0.5612;
        auto setup = makeSetup(50.0, 50.0, halfDepth);
        const double x0 = setup->lead->getRadLength();
        ShowerTally tally;
        tally.zFront = -halfDepth;
        tally.profile.assign(160, 0.0); // 160 x 0.25 X0 = 40 X0
        attachTally(setup->runManager, tally, x0);
        long tracks = 0;
        std::size_t maxStack = 0;
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < c.events; ++i) {
            Track primary(c.p, Vector3D(0, 0, -halfDepth - 0.5), Vector3D(0, 0, 1), c.energy);
            setup->runManager.simulateEvent(primary, 1.0 * Units::m, 10000000);
            tracks += setup->runManager.getNumberOfTracks();
            maxStack = std::max(maxStack, setup->runManager.getMaxStackSize());
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        double s0 = 0.0, s1 = 0.0, s2 = 0.0;
        for (std::size_t b = 0; b < tally.profile.size(); ++b) {
            const double t = (b + 0.5) * tally.binWidth;
            s0 += tally.profile[b];
            s1 += tally.profile[b] * t;
            s2 += tally.profile[b] * t * t;
        }
        const double mean = s1 / s0, var = s2 / s0 - mean * mean;
        const double b = mean / var, a = mean * b;
        const double tMaxPdg = std::log(c.energy / ec) + c.c;
        // peak: the largest 3-bin running mean, refined with a parabola through its neighbours
        std::size_t best = 1;
        auto smooth = [&](std::size_t i) { return tally.profile[i - 1] + tally.profile[i] + tally.profile[i + 1]; };
        for (std::size_t i = 2; i + 2 < tally.profile.size(); ++i) {
            if (smooth(i) > smooth(best)) best = i;
        }
        const double ym = smooth(best - 1), y0 = smooth(best), yp = smooth(best + 1);
        const double shift = 0.5 * (ym - yp) / (ym - 2.0 * y0 + yp);
        const double tPeak = (best + 0.5 + shift) * tally.binWidth;
        std::printf("%-13s %6d | %6.2f %6.3f | %7.2f %7.2f %7.2f | %8.2f%% %9.0f %9zu %7.2f\n", c.name, c.events, a, b,
                    (a - 1.0) / b, tPeak, tMaxPdg, 100.0 * s0 / (c.events * c.energy),
                    static_cast<double>(tracks) / c.events, maxStack, seconds / c.events);
        if (c.energy == 10.0 * Units::GeV && c.p == &electronDef) {
            shape10 = tally.profile;
            for (double& v : shape10) v /= c.events * c.energy * tally.binWidth;
        }
    }

    std::printf("\ne- 10 GeV: (1/E0) dE/dt per X0, MC vs the PDG gamma function (b = 0.5, PDG t_max):\n");
    const double b = 0.5, tMax = std::log(10.0 * Units::GeV / ec) - 0.5, a = b * tMax + 1.0;
    std::printf("%6s %10s %10s\n", "t", "MC", "PDG");
    for (int i = 2; i < 120 && !shape10.empty(); i += 8) {
        const double t = (i + 0.5) * 0.25;
        const double pdg = b * std::pow(b * t, a - 1.0) * std::exp(-b * t) / std::tgamma(a);
        std::printf("%6.2f %10.4f %10.4f\n", t, shape10[i], pdg);
    }
}

} // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "quick";
    RandomEngine::setSeed(12345);
    std::printf("[RandomEngine] Seed: %llu%s\n", static_cast<unsigned long long>(RandomEngine::getSeed()),
                quick ? "   (quick mode)" : "");

    const auto xcom = readXcom();
    const auto estar = readEstar();
    std::printf("Read %zu XCOM rows and %zu ESTAR rows from %s\n", xcom.size(), estar.size(), TOYSIM_DATA_DIR);

    testCrossSections(xcom);
    testAttenuation(quick ? 5000 : 200000);
    testSampling(quick ? 50000 : 400000);
    testKinematics(quick ? 20000 : 100000);
    testElectrons(estar);
    testEnergyConservation(quick ? 5 : 50);
    testShowerProfile(quick ? 50 : 2500);
    return 0;
}
