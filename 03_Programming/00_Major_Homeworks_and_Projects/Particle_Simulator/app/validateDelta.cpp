// ==============================================================================
// [STAGE 5] Validation of the delta rays.
//
//   TEST 1  The closed-form cross sections against a numerical integration of
//           dsigma/dT', and the restricted dE/dx against the textbook formulas
//           (restricted Bethe-Bloch, PDG eq. 34.7; restricted Berger-Seltzer for e-)
//   TEST 2  Sampled delta-ray spectra against dsigma/dT' (chi^2)
//   TEST 3  Kinematics of single collisions: energy and momentum conservation
//   TEST 4  Thin lead layer: the MEAN energy loss of a muon must not depend on
//           the cut, the number of delta rays must be sigma * L, and the energy
//           DEPOSITED in the layer drops (delta rays escape)
//   TEST 5  The mean range of a muon and the mean track length of an electron
//           must not depend on the cut
//   TEST 6  1 GeV electron showers: energy balance, track count, lateral spread
//
// Run:  ./ValidationDeltaApp          (all tests)
//       ./ValidationDeltaApp quick    (fewer events)
// ==============================================================================
#include "Annihilation.h"
#include "Box.h"
#include "Bremsstrahlung.h"
#include "ComptonScattering.h"
#include "Decay.h"
#include "GammaConversion.h"
#include "Ionization.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "ParticleTable.h"
#include "PhotoElectricEffect.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using PhysicalConstants::K;
using PhysicalConstants::m_e;

namespace {

const ParticleTable particleTable;

const Particle* find(const char* name) { return particleTable.findParticle(name); }

void header(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

std::shared_ptr<Material> makeLead() {
    DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

const std::shared_ptr<Material> lead = makeLead();

// Simpson in ln(T'):  integral f(T') dT' = integral f(T') T' d(ln T')
template <typename F>
double integrateLog(F f, double lo, double hi, int n = 4000) {
    const double a = std::log(lo), h = (std::log(hi) - a) / n;
    double s = 0.0;
    for (int i = 0; i <= n; ++i) {
        // exactly hi at the end: exp(log(hi)) can come out 1 ulp above Tmax, where dsigma = 0
        const double t = (i == n) ? hi : std::exp(a + i * h);
        s += ((i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0)) * f(t) * t;
    }
    return s * h / 3.0;
}

// Ionization with delta rays above `cut`
std::unique_ptr<Ionization> makeIonization(double cut) {
    auto ion = std::make_unique<Ionization>();
    ion->setFluctuationModel(FluctuationModel::landau);
    if (cut > 0.0) ion->setDeltaRayProduction(find("e-"), cut);
    ion->setVerboseLevel(0);
    return ion;
}

std::string cutLabel(double cut) {
    if (cut <= 0.0) return "off";
    char s[32];
    std::snprintf(s, sizeof(s), "%g MeV", cut / Units::MeV);
    return s;
}

// ==============================================================================
// TEST 1  Formulas
// ==============================================================================
// The closed forms of sigma and dE/dx_hard (Ionization.cpp) were integrated by
// hand; here the same integrals are done numerically from dsigma/dT' itself.
// Then the restricted dE/dx (unrestricted - hard) is compared with the formulas
// one finds in the literature, which were derived independently:
//   heavy: PDG eq. 34.7 (G4BetheBlochModel::ComputeDEDXPerVolume)
//     -dE/dx = (K/2) z^2 (Z/A) rho / beta^2 [ ln(2 m_e beta^2 gamma^2 T_up / I^2)
//              - beta^2 (1 + T_up/Tmax) - delta  (+ (T_up / 2E)^2, spin 1/2) ]
//   e-: restricted Berger-Seltzer (ICRU 37; G4MollerBhabhaModel), d = T_up / m_e:
//     -dE/dx = (K/2) (Z/A) rho / beta^2 [ ln(2 (tau+2) / (I/m_e)^2) - 1 - beta^2
//              + ln((tau-d) d) + tau/(tau-d) + (d^2/2 + (2 tau + 1) ln(1 - d/tau)) / gamma^2 - delta ]
double restrictedBetheBloch(const Ionization& ion, const Particle& p, double t, double cut) {
    const double mass = p.getMass(), e = t + mass, gamma = e / mass;
    const double bg2 = t * (t + 2.0 * mass) / (mass * mass), beta2 = bg2 / (gamma * gamma);
    const double tmax = ion.computeMaxSecondaryEnergy(p, t), tup = std::min(cut, tmax);
    const double i = lead->getMeanExcEn();
    double bracket = std::log(2.0 * m_e * bg2 * tup / (i * i)) - beta2 * (1.0 + tup / tmax)
                   - lead->DensityCorrection(0.5 * std::log10(bg2));
    if (p.getSpin() == 0.5) bracket += (tup / (2.0 * e)) * (tup / (2.0 * e));
    return 0.5 * K * p.getCharge() * p.getCharge() * lead->getZoverA() * lead->getDensity() / beta2 * bracket;
}

double restrictedBergerSeltzer(double t, double cut) {
    const double tau = t / m_e, gamma = tau + 1.0, bg2 = tau * (tau + 2.0), beta2 = bg2 / (gamma * gamma);
    const double i = lead->getMeanExcEn() / m_e, d = std::min(cut, 0.5 * t) / m_e;
    const double bracket = std::log(2.0 * (tau + 2.0) / (i * i)) - 1.0 - beta2 + std::log((tau - d) * d) + tau / (tau - d)
                         + (0.5 * d * d + (2.0 * tau + 1.0) * std::log(1.0 - d / tau)) / (gamma * gamma)
                         - lead->DensityCorrection(0.5 * std::log10(bg2));
    return 0.5 * K * lead->getZoverA() * lead->getDensity() / beta2 * bracket;
}

void testFormulas() {
    header("TEST 1  sigma and dE/dx_hard: closed form vs numerical integral (lead)");
    std::printf("%-6s %9s %8s | %10s %9s | %10s %9s | %8s | %9s\n", "", "T [MeV]", "cut", "sigma[1/cm]", "rel.diff",
                "hard[MeV/cm]", "rel.diff", "restr/all", "vs book");
    const Ionization ion; // the formulas do not need delta rays switched on
    struct Case { const char* name; double t; double cut; };
    for (const Case& c : {Case{"mu+", 100.0, 1.0}, Case{"mu+", 1000.0, 1.0}, Case{"mu+", 1.0e5, 1.0}, Case{"mu+", 1000.0, 0.1},
                          Case{"pi+", 1000.0, 1.0}, Case{"tau-", 1.0e4, 1.0}, Case{"e-", 3.0, 1.0}, Case{"e-", 100.0, 1.0},
                          Case{"e-", 1.0e4, 0.1}, Case{"e+", 3.0, 1.0}, Case{"e+", 100.0, 1.0}, Case{"e+", 1.0e4, 0.1}}) {
        const Particle& p = *find(c.name);
        const double t = c.t * Units::MeV, cut = c.cut * Units::MeV;
        const double tmax = std::min(ion.computeMaxSecondaryEnergy(p, t), t);
        auto dsig = [&](double tp) { return ion.computeDeltaRayDifferentialCrossSection(p, t, tp, *lead); };
        const double sigma = ion.computeDeltaRayCrossSection(p, t, *lead, cut);
        const double hard = ion.computeDeltaRayEnergyLoss(p, t, *lead, cut);
        const double sigmaNum = integrateLog(dsig, cut, tmax);
        const double hardNum = integrateLog([&](double tp) { return tp * dsig(tp); }, cut, tmax);
        const double restricted = ion.computeRestrictedDEDX(p, t, *lead, cut);
        const double book = (std::abs(p.getPdgId()) == 11) ? (p.getPdgId() == 11 ? restrictedBergerSeltzer(t, cut) : -1.0)
                                                           : restrictedBetheBloch(ion, p, t, cut);
        char bookText[16] = "      -";
        if (book > 0.0) std::snprintf(bookText, sizeof(bookText), "%9.1e", restricted / book - 1.0);
        std::printf("%-6s %9g %8s | %10.4g %9.1e | %10.4g %9.1e | %8.4f | %s\n", c.name, c.t, cutLabel(cut).c_str(),
                    sigma * Units::cm, sigma / sigmaNum - 1.0, hard / (Units::MeV / Units::cm), hard / hardNum - 1.0,
                    restricted / ion.computeDEDX(p, t, *lead), bookText);
    }
    std::printf("rel.diff: closed form / Simpson - 1 (rounding, <= ~1e-12). restr/all: share of dE/dx that\n");
    std::printf("stays continuous. vs book: restricted dE/dx / literature formula - 1 (e+: no closed form quoted).\n");
}

// ==============================================================================
// TEST 2  Sampled spectra
// ==============================================================================
void testSpectra(int nSamples) {
    header("TEST 2  Sampled delta-ray spectra vs dsigma/dT' (40 bins in ln T')");
    Ionization ion;
    struct Case { const char* name; double t; double cut; };
    for (const Case& c : {Case{"mu+", 1.0e4, 1.0}, Case{"mu+", 300.0, 1.0}, Case{"e-", 100.0, 1.0}, Case{"e-", 3.0, 1.0},
                          Case{"e+", 100.0, 1.0}, Case{"e+", 3.0, 1.0}}) {
        const Particle& p = *find(c.name);
        const double t = c.t * Units::MeV, cut = c.cut * Units::MeV;
        const double tmax = std::min(ion.computeMaxSecondaryEnergy(p, t), t);
        const int nBins = 40;
        const double lo = std::log(cut), hi = std::log(tmax);
        std::vector<double> observed(nBins, 0.0), expected(nBins, 0.0);
        double sum = 0.0, sum2 = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            const double tp = ion.sampleDeltaRayEnergy(p, t, cut);
            sum += tp;
            sum2 += tp * tp;
            const int b = static_cast<int>((std::log(tp) - lo) / (hi - lo) * nBins);
            observed[std::min(std::max(b, 0), nBins - 1)] += 1.0;
        }
        auto dsig = [&](double tp) { return ion.computeDeltaRayDifferentialCrossSection(p, t, tp, *lead); };
        const double total = ion.computeDeltaRayCrossSection(p, t, *lead, cut);
        for (int b = 0; b < nBins; ++b) {
            expected[b] = nSamples * integrateLog(dsig, std::exp(lo + (hi - lo) * b / nBins), std::exp(lo + (hi - lo) * (b + 1) / nBins), 200) / total;
        }
        double chi2 = 0.0;
        int ndf = 0;
        for (int b = 0; b < nBins; ++b) {
            if (expected[b] < 5.0) continue;
            chi2 += (observed[b] - expected[b]) * (observed[b] - expected[b]) / expected[b];
            ++ndf;
        }
        const double meanTheory = ion.computeDeltaRayEnergyLoss(p, t, *lead, cut) / total;
        const double mean = sum / nSamples, error = std::sqrt((sum2 / nSamples - mean * mean) / nSamples);
        std::printf("  %-4s T = %6g MeV, T' in [%g, %.4g] MeV: chi2/ndf = %5.1f / %2d = %4.2f | <T'> %.4f vs %.4f MeV (%+.2f)\n",
                    c.name, c.t, c.cut, tmax / Units::MeV, chi2, ndf - 1, chi2 / std::max(ndf - 1, 1),
                    mean / Units::MeV, meanTheory / Units::MeV, (mean - meanTheory) / error);
    }
    std::printf("  ndf = bins - 1 (the total is fixed). (...): pull of <T'> = mean energy of a delta ray.\n");
}

// ==============================================================================
// TEST 3  Kinematics of single collisions
// ==============================================================================
void testKinematics(int n) {
    header("TEST 3  Single collisions: T - T_after - T' and |p - p_after - p'| / p (max over events)");
    auto ion = makeIonization(1.0 * Units::MeV);
    struct Case { const char* name; double t; };
    for (const Case& c : {Case{"mu+", 1000.0}, Case{"mu-", 1.0e5}, Case{"tau-", 1.0e4}, Case{"e-", 10.0}, Case{"e+", 3.0},
                          Case{"e+", 1000.0}}) {
        const Particle* p = find(c.name);
        const double t = c.t * Units::MeV;
        double maxDe = 0.0, maxDp = 0.0, sumCos = 0.0, sumTp = 0.0;
        int stopped = 0;
        for (int i = 0; i < n; ++i) {
            Track track(p, Vector3D(), Vector3D(0.6, 0.0, 0.8), t);
            const Vector3D p0 = track.getMomentum();
            Step step;
            std::vector<Track> secondaries;
            ion->postStepDoIt(track, *lead, step, secondaries);
            const Track& d = secondaries.at(0);
            maxDe = std::max(maxDe, std::abs(t - track.getKineticEnergy() - d.getKineticEnergy()) / t);
            maxDp = std::max(maxDp, (p0 - track.getMomentum() - d.getMomentum()).mag() / p0.mag());
            sumCos += d.getDirection().dot(p0.unit());
            sumTp += d.getKineticEnergy();
            if (track.getKineticEnergy() <= 0.0) ++stopped;
        }
        std::printf("  %-4s %8g MeV | max dE/T %8.1e | max dp/p %8.1e | <T'> %7.3f MeV  <cos theta'> %.4f | stopped %d\n",
                    c.name, c.t, maxDe, maxDp, sumTp / n / Units::MeV, sumCos / n, stopped);
    }
    std::printf("  [INFO] cos theta' = T'(E + m_e)/(p' p): a 1 MeV delta ray leaves at ~ 34 deg whatever the\n");
    std::printf("         projectile; only T' close to Tmax goes forward. e+ 3 MeV: eps = 1 gives all of T away.\n");
}

// ==============================================================================
// Full simulation: a lead box in a vacuum world 1 cm larger on each side
// ==============================================================================
struct Setup {
    RunManager runManager;
    Ionization* ionization = nullptr; // owned by the RunManager
};
std::unique_ptr<Setup> makeSetup(double hx, double hy, double hz, double deltaCut) {
    auto setup = std::make_unique<Setup>();
    auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
    RunManager& rm = setup->runManager;
    rm.setWorldVolume(std::make_unique<Box>("World", hx + 1.0, hy + 1.0, hz + 1.0, vacuum));
    rm.addVolume(std::make_unique<Box>("Lead", hx, hy, hz, lead));
    const Particle* gamma = find("gamma");
    const Particle* electron = find("e-");
    rm.addPhysicsProcess(std::make_unique<MultipleScattering>());
    auto ionization = makeIonization(deltaCut);
    setup->ionization = ionization.get();
    rm.addPhysicsProcess(std::move(ionization));
    rm.addPhysicsProcess(std::make_unique<Decay>());
    rm.addPhysicsProcess(std::make_unique<Bremsstrahlung>(gamma));
    rm.addPhysicsProcess(std::make_unique<Annihilation>(gamma));
    rm.addPhysicsProcess(std::make_unique<PhotoElectricEffect>(electron));
    rm.addPhysicsProcess(std::make_unique<ComptonScattering>(electron));
    rm.addPhysicsProcess(std::make_unique<GammaConversion>(electron, find("e+")));
    rm.setEnergyCut(100.0 * Units::keV);
    rm.setVerboseLevel(0);
    return setup;
}

// Per-event bookkeeping. The primary is the first track the RunManager follows
// (it tracks one particle to its end before it takes the next from the stack).
struct EventRecord {
    bool primaryActive = true;
    double primaryExitEnergy = -1.0;  // kinetic energy when the primary left the world
    double primaryLength = 0.0;       // path length of the primary INSIDE THE LEAD (a backscattered
                                      // electron can fly 30 cm through the vacuum gap at a grazing angle)
    int primaryDeltas = 0;
    int deltas = 0;
    double deposited = 0.0, escaped = 0.0, escapedOther = 0.0;
    double sumRE = 0.0;               // integral r dE (lateral spread)
    int nConv = 0, nAnnih = 0;
    void clear() { *this = EventRecord(); }
};
void attach(Setup& setup, EventRecord& ev, const Vector3D& leadHalf) {
    setup.runManager.setUserSteppingAction([&ev, leadHalf](const Track& track, const Step& step) {
        const Vector3D& pre = step.preStepPosition;
        const bool inLead = std::abs(pre.x) <= leadHalf.x && std::abs(pre.y) <= leadHalf.y && std::abs(pre.z) <= leadHalf.z;
        ev.deposited += step.totalEnergyDeposit;
        const Vector3D mid = (step.preStepPosition + step.postStepPosition) * 0.5;
        ev.sumRE += std::sqrt(mid.x * mid.x + mid.y * mid.y) * step.totalEnergyDeposit;
        if (step.processDefinedStep == "conv") ++ev.nConv;
        if (step.processDefinedStep == "annihil") ++ev.nAnnih;
        const int nDelta = (step.processDefinedStep == "Ionization" && step.secondaryInCurrentStep != nullptr)
                         ? static_cast<int>(step.secondaryInCurrentStep->size()) : 0;
        ev.deltas += nDelta;
        if (step.leftWorld) {
            ev.escaped += track.getKineticEnergy();
            if (!ev.primaryActive) ev.escapedOther += track.getKineticEnergy();
        }
        if (ev.primaryActive) {
            ev.primaryDeltas += nDelta;
            if (inLead) ev.primaryLength += step.stepLength;
            if (step.leftWorld) ev.primaryExitEnergy = track.getKineticEnergy();
            if (track.getTrackStatus() == TrackStatus::trackKill) ev.primaryActive = false;
        }
    });
}

struct Mean {
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    void add(double x) { sum += x; sum2 += x * x; ++n; }
    double mean() const { return sum / n; }
    double error() const { return std::sqrt(std::max(0.0, sum2 / n - mean() * mean()) / n); }
};

const std::vector<double> cuts = {0.0, 10.0 * Units::MeV, 1.0 * Units::MeV, 0.1 * Units::MeV};

// ==============================================================================
// TEST 4  Thin layer
// ==============================================================================
// 1 GeV mu+ through 1 cm of lead. The muon's own mean loss T0 - T_exit is the
// unrestricted dE/dx * L whatever the cut. The energy DEPOSITED in the layer is
// not: a delta ray of up to Tmax = 99 MeV can carry its energy out of a 1 cm layer.
void testThinLayer(int nEvents) {
    header("TEST 4  mu+ 1 GeV through 1 cm of lead: mean loss vs deposited energy");
    const Particle* mu = find("mu+");
    const double t0 = 1.0 * Units::GeV, thickness = 1.0 * Units::cm;
    const Ionization reference;
    const double expectedLoss = reference.computeDEDX(*mu, t0, *lead) * thickness;
    std::printf("%-8s | %15s %8s | %7s %7s | %9s %9s | %9s %5s\n", "cut", "<T0 - T_exit>", "pull", "N_delta", "sigma L",
                "deposited", "esc.other", "max imbal.", "viol.");
    for (double cut : cuts) {
        auto setup = makeSetup(25.0, 25.0, 0.5 * thickness, cut);
        EventRecord ev;
        attach(*setup, ev, Vector3D(25.0, 25.0, 0.5 * thickness));
        Mean loss, deltas, dep, other;
        double maxImbalance = 0.0;
        for (int i = 0; i < nEvents; ++i) {
            ev.clear();
            setup->runManager.simulateEvent(Track(mu, Vector3D(0, 0, -1.0), Vector3D(0, 0, 1), t0), 1.0 * Units::m, 100000);
            loss.add(t0 - ev.primaryExitEnergy);
            deltas.add(ev.primaryDeltas);
            dep.add(ev.deposited);
            other.add(ev.escapedOther);
            maxImbalance = std::max(maxImbalance, std::abs(ev.deposited + ev.escaped - t0 - 2.0 * m_e * (ev.nAnnih - ev.nConv)) / t0);
        }
        const double sigmaL = cut > 0.0 ? reference.computeDeltaRayCrossSection(*mu, t0, *lead, cut) * thickness : 0.0;
        std::printf("%-8s | %8.3f+-%.3f %+8.2f | %7.3f %7.3f | %9.3f %9.3f | %9.1e %5ld\n", cutLabel(cut).c_str(), loss.mean() / Units::MeV,
                    loss.error() / Units::MeV, (loss.mean() - expectedLoss) / loss.error(), deltas.mean(), sigmaL,
                    dep.mean() / Units::MeV, other.mean() / Units::MeV, maxImbalance, setup->ionization->getSigmaMaxViolations());
    }
    std::printf("Energies in MeV. Expected <T0 - T_exit> = dE/dx(1 GeV) * 1 cm = %.3f MeV (pull vs this; msc makes\n", expectedLoss / Units::MeV);
    std::printf("the path ~0.01 %% longer). sigma L: expected delta rays above the cut. esc.other: energy that\n");
    std::printf("leaves the layer with particles other than the muon (delta rays and their photons).\n");
    std::printf("viol.: steps where sigma at the end of the step was above the sampled sigmaMax (should be ~0).\n");
}

// ==============================================================================
// TEST 5  Mean range / track length
// ==============================================================================
void testRange(int nEvents) {
    header("TEST 5  Mean path length of the primary in the lead vs the cut (block 60 x 60 x 60 cm)");
    std::printf("%-12s %-8s | %16s %8s | %8s %8s %8s | %5s\n", "primary", "cut", "<track length>", "vs off", "N_delta", "expected",
                "pull", "viol.");
    struct Case { const char* name; double t; const char* label; };
    for (const Case& c : {Case{"mu+", 200.0 * Units::MeV, "mu+ 200 MeV"}, Case{"e-", 20.0 * Units::MeV, "e- 20 MeV"}}) {
        const Particle* p = find(c.name);
        Mean off;
        for (double cut : cuts) {
            auto setup = makeSetup(30.0, 30.0, 30.0, cut);
            EventRecord ev;
            attach(*setup, ev, Vector3D(30.0, 30.0, 30.0));
            Mean length, deltas;
            for (int i = 0; i < nEvents; ++i) {
                ev.clear();
                setup->runManager.simulateEvent(Track(p, Vector3D(0, 0, -30.5), Vector3D(0, 0, 1), c.t), 1.0 * Units::m, 1000000);
                length.add(ev.primaryLength);
                deltas.add(ev.primaryDeltas);
            }
            if (cut <= 0.0) off = length;
            const double diff = length.mean() - off.mean();
            const double err = std::sqrt(length.error() * length.error() + off.error() * off.error());
            char pull[16] = "       -";
            if (cut > 0.0) std::snprintf(pull, sizeof(pull), "%+8.2f", diff / err);
            // Expected number of delta rays of a muon (no radiative loss, so dx = dT / S on average):
            //   N = integral sigma(T) dx = integral_0^T0 sigma(T) / S(T) dT
            char expected[40] = "       -        -";
            if (std::abs(p->getPdgId()) == 13 && cut > 0.0 && deltas.mean() > 0.0) {
                const Ionization formulas;
                const double n = integrateLog([&](double t) { return formulas.computeDeltaRayCrossSection(*p, t, *lead, cut)
                                                                     / formulas.computeDEDX(*p, t, *lead); }, 1.0 * Units::keV, c.t, 20000);
                std::snprintf(expected, sizeof(expected), "%8.3f %+8.2f", n, (deltas.mean() - n) / deltas.error());
            }
            std::printf("%-12s %-8s | %9.4f+-%.4f %s | %8.3f %s | %5ld\n", c.label, cutLabel(cut).c_str(), length.mean(), length.error(), pull,
                        deltas.mean(), expected, setup->ionization->getSigmaMaxViolations());
        }
        if (std::abs(p->getPdgId()) == 13) {
            Ionization csda;
            csda.setVerboseLevel(0);
            std::printf("%-12s CSDA range from the unrestricted table: %.4f cm\n", c.label, csda.getRange(*p, c.t, *lead));
        }
    }
    std::printf("Path lengths in the lead, cm. vs off: (cut - off) / error. Electrons: until the 100 keV energy cut;\n");
    std::printf("in Moller scattering the faster electron is kept as the primary (Tmax = T/2). expected: the muon's\n");
    std::printf("mean number of delta rays, integral sigma(T)/S(T) dT (S: unrestricted dE/dx).\n");
}

// ==============================================================================
// TEST 6  Showers
// ==============================================================================
void testShower(int nEvents) {
    header("TEST 6  e- 1 GeV showers in a lead block 20 x 20 x 20 cm");
    std::printf("%-8s | %9s %10s | %8s %8s | %8s %7s | %5s\n", "cut", "deposited", "max imbal.", "tracks", "N_delta", "<r> [cm]",
                "s/evt", "viol.");
    const Particle* e = find("e-");
    const double t0 = 1.0 * Units::GeV;
    for (double cut : cuts) {
        auto setup = makeSetup(10.0, 10.0, 10.0, cut);
        EventRecord ev;
        attach(*setup, ev, Vector3D(10.0, 10.0, 10.0));
        Mean dep, tracks, deltas, radius;
        double maxImbalance = 0.0;
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < nEvents; ++i) {
            ev.clear();
            setup->runManager.simulateEvent(Track(e, Vector3D(0, 0, -10.5), Vector3D(0, 0, 1), t0), 1.0 * Units::m, 100000);
            dep.add(ev.deposited / t0);
            tracks.add(setup->runManager.getNumberOfTracks());
            deltas.add(ev.deltas);
            radius.add(ev.sumRE / ev.deposited);
            maxImbalance = std::max(maxImbalance, std::abs(ev.deposited + ev.escaped - t0 - 2.0 * m_e * (ev.nAnnih - ev.nConv)) / t0);
        }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::printf("%-8s | %8.2f%% %10.1e | %8.0f %8.1f | %8.3f %7.3f | %5ld\n", cutLabel(cut).c_str(), 100.0 * dep.mean(), maxImbalance,
                    tracks.mean(), deltas.mean(), radius.mean(), seconds / nEvents, setup->ionization->getSigmaMaxViolations());
    }
    std::printf("<r>: energy-weighted mean distance of the deposits from the shower axis.\n");
}

} // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "quick";
    RandomEngine::setSeed(12345);
    std::printf("[RandomEngine] Seed: %llu%s\n", static_cast<unsigned long long>(RandomEngine::getSeed()),
                quick ? "  (quick mode)" : "");

    const int f = quick ? 10 : 1;
    testFormulas();
    testSpectra(1000000 / f);
    testKinematics(100000 / f);
    testThinLayer(20000 / f);
    testRange(4000 / f);
    testShower(300 / f);
    return 0;
}
