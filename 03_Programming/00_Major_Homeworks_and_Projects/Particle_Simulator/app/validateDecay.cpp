// ==============================================================================
// Validation of the decays.
//
//   TEST 1  Muon decay at rest: Michel spectrum of the e+, the two neutrino spectra
//   TEST 2  Leptonic tau decays: Michel spectrum with the lepton mass (tau -> mu nu nu)
//   TEST 3  Branching ratios: how often each channel is chosen
//   TEST 4  Kinematics of single decays: two-body momenta, energy and momentum
//           conservation after the boost, invariant mass of the products
//   TEST 5  Decay law in vacuum: proper time and decay length vs tau, beta gamma c tau
//   TEST 6  Decay law with energy loss: a pi+ / mu+ that slows down and STOPS in lead
//           still decays at the proper time drawn at its start
//   TEST 7  Energy conservation in full events, neutrinos counted as escaping
//
// Run:  ./ValidationDecayApp          (all tests, ~1 min)
//       ./ValidationDecayApp quick    (fewer events, a few seconds)
// ==============================================================================
#include "Annihilation.h"
#include "Box.h"
#include "Bremsstrahlung.h"
#include "ComptonScattering.h"
#include "Decay.h"
#include "DecayTable.h"
#include "GammaConversion.h"
#include "Ionization.h"
#include "LeptonicDecayChannel.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "ParticleTable.h"
#include "PhotoElectricEffect.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using PhysicalConstants::m_e;

namespace {

const ParticleTable particleTable;

const Particle* find(const char* name) { return particleTable.findParticle(name); }

bool isNeutrino(const Particle& p) {
    const int pdg = std::abs(p.getPdgId());
    return pdg == 12 || pdg == 14 || pdg == 16;
}

void header(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
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
    std::printf("  %-44s chi2/ndf = %7.1f / %3d = %5.2f\n", label, chi2, ndf, chi2 / std::max(ndf, 1));
}

// N * (integral of the normalised pdf over each bin), Simpson's rule
template <typename Pdf>
std::vector<double> expectedCounts(Pdf pdf, double lo, double hi, int nBins, int nSamples) {
    std::vector<double> integral(nBins, 0.0);
    double total = 0.0;
    for (int b = 0; b < nBins; ++b) {
        const double a = lo + (hi - lo) * b / nBins;
        const double c = lo + (hi - lo) * (b + 1) / nBins;
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

void fill(std::vector<double>& h, double x, double lo, double hi) {
    const int b = static_cast<int>((x - lo) / (hi - lo) * h.size());
    if (b >= 0 && b < static_cast<int>(h.size())) h[b] += 1.0;
}

std::shared_ptr<Material> makeLead() {
    DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

// A lead box (half-sizes hx, hy, hz) in a vacuum world 1 cm larger on every side,
// the full physics list of stages 2-4
struct Setup {
    RunManager runManager;
    std::shared_ptr<Material> lead;
};
std::unique_ptr<Setup> makeSetup(double hx, double hy, double hz) {
    auto setup = std::make_unique<Setup>();
    setup->lead = makeLead();
    auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
    RunManager& rm = setup->runManager;
    rm.setWorldVolume(std::make_unique<Box>("World", hx + 1.0, hy + 1.0, hz + 1.0, vacuum));
    rm.addVolume(std::make_unique<Box>("Lead", hx, hy, hz, setup->lead));
    const Particle* gamma = find("gamma");
    const Particle* electron = find("e-");
    rm.addPhysicsProcess(std::make_unique<MultipleScattering>());
    rm.addPhysicsProcess(std::make_unique<Ionization>());
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

// The daughters' rest-frame four-momenta of one decay of `parent`, with the
// channel that was chosen
struct RestFrameDecay {
    const DecayChannel* channel;
    std::vector<LorentzVector> products;
};
RestFrameDecay decayAtRest(const Particle* parent) {
    const DecayChannel* channel = parent->getDecayTable()->selectADecayChannel();
    return {channel, channel->decayIt(parent->getMass())};
}

// ==============================================================================
// TEST 1  Muon decay at rest
// ==============================================================================
// mu+ -> e+ nu_e anti_nu_mu. x = E / E_max.
//   e+:          Michel, sqrt(x^2 - x0^2)(3x - 2x^2 - x0^2) ~ x^2 (3 - 2x),  <x> = 0.7
//   anti_nu_mu:  same shape as the e+ (it is the one paired with the e+ in |M|^2)
//   nu_e:        12 x^2 (1 - x),                                            <x> = 0.6
// Mean energies, in units of m_mu: 0.35 + 0.30 + 0.35 = 1.
void testMuonDecay(int n) {
    header("TEST 1  mu+ -> e+ nu_e anti_nu_mu at rest: spectra (x = E / E_max)");
    const Particle* mu = find("mu+");
    const double M = mu->getMass();
    const double mE = find("e+")->getMass();
    const double eMaxLepton = (M * M + mE * mE) / (2.0 * M);
    const double eMaxNu = (M * M - mE * mE) / (2.0 * M);
    const double x0 = mE / eMaxLepton;
    const int nBins = 50;
    std::vector<double> hE(nBins, 0.0), hNu2(nBins, 0.0), hNu3(nBins, 0.0);
    double sumE[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i) {
        const RestFrameDecay d = decayAtRest(mu);
        fill(hE, d.products[0].E() / eMaxLepton, 0.0, 1.0);
        fill(hNu2, d.products[1].E() / eMaxNu, 0.0, 1.0);
        fill(hNu3, d.products[2].E() / eMaxNu, 0.0, 1.0);
        for (int k = 0; k < 3; ++k) sumE[k] += d.products[k].E();
    }
    printChi2("e+ vs Michel (rho = 3/4, with m_e)", hE,
              expectedCounts([x0](double x) { return LeptonicDecayChannel::michelSpectrum(x, x0); }, 0.0, 1.0, nBins, n));
    printChi2("anti_nu_mu vs x^2 (3 - 2x)", hNu3, expectedCounts([](double x) { return x * x * (3.0 - 2.0 * x); }, 0.0, 1.0, nBins, n));
    printChi2("nu_e vs x^2 (1 - x)", hNu2, expectedCounts([](double x) { return x * x * (1.0 - x); }, 0.0, 1.0, nBins, n));
    std::printf("  mean energy / m_mu:   e+ %.4f (0.35)   nu_e %.4f (0.30)   anti_nu_mu %.4f (0.35)   sum %.6f\n",
                sumE[0] / n / M, sumE[1] / n / M, sumE[2] / n / M, (sumE[0] + sumE[1] + sumE[2]) / n / M);
    std::printf("  -> the e+ carries 35 %% of the muon's rest energy on average, the neutrinos 65 %%\n");
    std::printf("     (the old toy model gave the electron 80 %% of the KINETIC energy, 0 at rest).\n");
    std::printf("  Michel spectrum, normalised dN/dx and MC:\n     x   formula      MC\n");
    const auto expE = expectedCounts([x0](double x) { return LeptonicDecayChannel::michelSpectrum(x, x0); }, 0.0, 1.0, nBins, n);
    for (int b = 2; b < nBins; b += 6) {
        std::printf("  %5.2f %8.4f %8.4f\n", (b + 0.5) / nBins, expE[b] * nBins / n, hE[b] * nBins / n);
    }
}

// ==============================================================================
// TEST 2  Leptonic tau decays
// ==============================================================================
// Same matrix element, but for tau -> mu nu nu the lepton mass is not small:
// x0 = m_mu / E_max = 0.12. The Michel formula with x0 tests the mass terms.
void testTauLeptonic(int n) {
    header("TEST 2  tau- -> l- anti_nu_l nu_tau at rest: lepton spectrum vs Michel with mass");
    const Particle* tau = find("tau-");
    const double M = tau->getMass();
    const DecayTable& table = *tau->getDecayTable();
    for (int c = 0; c < 2; ++c) {
        const DecayChannel* channel = table.getDecayChannel(c);
        const double m1 = channel->getDaughter(0)->getMass();
        const double eMax = (M * M + m1 * m1) / (2.0 * M);
        const double x0 = m1 / eMax;
        const int nBins = 50;
        std::vector<double> h(nBins, 0.0);
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const std::vector<LorentzVector> p = channel->decayIt(M);
            fill(h, p[0].E() / eMax, 0.0, 1.0);
            sum += p[0].E();
        }
        char label[80];
        std::snprintf(label, sizeof(label), "%s (x0 = %.4f)", channel->getDescription().c_str(), x0);
        printChi2(label, h, expectedCounts([x0](double x) { return LeptonicDecayChannel::michelSpectrum(x, x0); }, 0.0, 1.0, nBins, n));
        std::printf("  %-44s <E_lepton> / m_tau = %.4f\n", "", sum / n / M);
    }
}

// ==============================================================================
// TEST 3  Branching ratios
// ==============================================================================
void testBranchingRatios(int n) {
    header("TEST 3  Branching ratios: fraction of decays in each channel");
    std::printf("  %-36s %9s %9s %8s\n", "channel", "BR", "MC", "pull");
    for (const char* name : {"tau-", "pi+"}) {
        const DecayTable& table = *find(name)->getDecayTable();
        std::vector<int> count(table.entries(), 0);
        for (int i = 0; i < n; ++i) {
            const DecayChannel* c = table.selectADecayChannel();
            for (int k = 0; k < table.entries(); ++k) {
                if (table.getDecayChannel(k) == c) ++count[k];
            }
        }
        double hadronicBR = 0.0, hadronicMC = 0.0;
        for (int k = 0; k < table.entries(); ++k) {
            const DecayChannel* c = table.getDecayChannel(k);
            const double br = c->getBR() / table.getSumOfBR();
            const double mc = static_cast<double>(count[k]) / n;
            const double sigma = std::sqrt(br * (1.0 - br) / n);
            std::printf("  %-36s %9.5f %9.5f %+8.2f\n", c->getDescription().c_str(), br, mc, (mc - br) / sigma);
            if (std::strcmp(name, "tau-") == 0 && k >= 2) {
                hadronicBR += br;
                hadronicMC += mc;
            }
        }
        if (std::strcmp(name, "tau-") == 0) {
            std::printf("  %-36s %9.5f %9.5f\n", "tau- hadronic (sum)", hadronicBR, hadronicMC);
        }
    }
    std::printf("  pull = (MC - BR) / binomial sigma: should be within ~ +-2 (+-3 for a few of them).\n");
}

// ==============================================================================
// TEST 4  Kinematics of single decays
// ==============================================================================
// The Decay process itself is called on a Track (in flight, so the boost is
// tested too). The daughters must have, together, the parent's energy and
// momentum, and their invariant mass must be the parent mass.
void testKinematics(int n) {
    header("TEST 4  Kinematics of single decays");
    // Two-body momenta at rest: fixed by energy-momentum conservation
    const Particle* piPlus = find("pi+");
    const Particle* piZero = find("pi0");
    const double pMu = DecayChannel::twoBodyMomentum(piPlus->getMass(), find("mu+")->getMass(), 0.0);
    std::printf("  pi+ -> mu+ nu_mu at rest: p = %.4f MeV (PDG 29.79), T_mu = %.4f MeV (4.12)\n",
                pMu / Units::MeV, (std::sqrt(pMu * pMu + std::pow(find("mu+")->getMass(), 2)) - find("mu+")->getMass()) / Units::MeV);
    double maxDev = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const RestFrameDecay d = decayAtRest(piZero);
        maxDev = std::max(maxDev, std::abs(d.products[0].E() - 0.5 * piZero->getMass()) / Units::MeV);
    }
    std::printf("  pi0 -> gamma gamma at rest: E_gamma = %.4f MeV (m_pi0 / 2), max deviation %.1e MeV\n\n",
                0.5 * piZero->getMass() / Units::MeV, maxDev);

    Decay decay;
    decay.setVerboseLevel(0);
    Material vacuum("Vacuum", 0.0);
    std::printf("  %-8s %10s | %12s %12s %12s\n", "parent", "p [GeV]", "max |dE|/E", "max |dp|/E", "max |dm|/m");
    for (const char* name : {"mu-", "mu+", "pi+", "pi-", "pi0", "tau-", "tau+"}) {
        const Particle* parent = find(name);
        for (double pGeV : {0.0, 0.1, 10.0, 1000.0}) {
            double maxE = 0.0, maxP = 0.0, maxM = 0.0;
            const double m = parent->getMass();
            const double kinetic = std::sqrt(pGeV * pGeV + m * m) - m;
            for (int i = 0; i < n; ++i) {
                // a random direction, so that the boost is not always along z
                const double ct = 2.0 * RandomEngine::uniNumGen() - 1.0, phi = 2.0 * M_PI * RandomEngine::uniNumGen();
                const Vector3D dir(std::sqrt(1.0 - ct * ct) * std::cos(phi), std::sqrt(1.0 - ct * ct) * std::sin(phi), ct);
                Track track(parent, Vector3D(), dir, kinetic);
                Step step;
                std::vector<Track> secondaries;
                decay.postStepDoIt(track, vacuum, step, secondaries);
                LorentzVector sum;
                for (const Track& s : secondaries) sum = sum + LorentzVector(s.getMomentum(), s.getTotalEnergy());
                const double e0 = track.getTotalEnergy();
                maxE = std::max(maxE, std::abs(sum.E() - e0) / e0);
                maxP = std::max(maxP, (sum.vect() - track.getMomentum()).mag() / e0);
                // invariant mass from E^2 - p^2 loses digits when E >> m: (E - p)(E + p) instead
                const double pSum = sum.vect().mag();
                const double mSum = std::sqrt(std::max(0.0, (sum.E() - pSum) * (sum.E() + pSum)));
                maxM = std::max(maxM, std::abs(mSum - m) / m);
            }
            std::printf("  %-8s %10.1f | %12.1e %12.1e %12.1e\n", name, pGeV, maxE, maxP, maxM);
        }
    }
    std::printf("  Only rounding (~1e-16 for E and p). The mass from E^2 - p^2 has a larger rounding\n");
    std::printf("  error when E >> m (p = 1000 GeV, gamma up to 1e4): ~1e-16 * gamma^2.\n");
}

// ==============================================================================
// TEST 5  Decay law in vacuum
// ==============================================================================
// The primary's proper time at the decay must be exponential with mean tau,
// and (no energy loss in vacuum) its flight length exponential with mean
// beta gamma c tau. The primary is the first track of the event and nothing
// else exists before it decays, so the first "Decay" step of an event is its decay.
void testDecayLawVacuum(int n) {
    header("TEST 5  Decay law in vacuum: proper time and decay length of the primary");
    std::printf("  %-6s %8s | %10s %10s %7s | %12s %12s %7s\n", "", "p [GeV]", "<tau> [ns]", "PDG", "pull",
                "<L> [cm]", "bgc tau", "pull");
    struct Case { const char* name; double pGeV; };
    // mu+ at 10 MeV/c: decay length 62 m. At 1 GeV/c it would be 6 km, and the
    // world would have to be ~100 km (see the kCarTolerance note below).
    for (const Case& c : {Case{"pi+", 1.0}, Case{"mu+", 0.01}, Case{"tau-", 10.0}, Case{"pi0", 1.0}}) {
        const Particle* parent = find(c.name);
        const double tau = parent->getLifetime();
        const double m = parent->getMass();
        const double kinetic = std::sqrt(c.pGeV * c.pGeV + m * m) - m;
        const double meanLength = c.pGeV / m * PhysicalConstants::c * tau;

        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        // 40 mean decay lengths: the primary leaves before decaying with probability
        // e^-40. Not more: kCarTolerance = 1e-10 cm is an ABSOLUTE push across a
        // boundary, and at |x| > ~4e5 cm (4 km) it is below the rounding of x
        // (4e5 * 2.2e-16 = 1e-10): x + kCarTolerance == x, and a particle on the
        // world wall would never get out (Geant4 has the same limit for huge worlds).
        const double half = 40.0 * meanLength;
        rm.setWorldVolume(std::make_unique<Box>("World", half, half, half, vacuum));
        rm.addPhysicsProcess(std::make_unique<Decay>());
        rm.setVerboseLevel(0);

        bool primaryDone = false;
        double sumTau = 0.0, sumL = 0.0;
        int decays = 0;
        const int nBins = 40;
        std::vector<double> h(nBins, 0.0);
        rm.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (primaryDone || step.processDefinedStep != "Decay") return;
            primaryDone = true;
            ++decays;
            sumTau += track.getProperTime();
            sumL += track.getTrackLength();
            fill(h, track.getProperTime() / tau, 0.0, 5.0);
        });
        for (int i = 0; i < n; ++i) {
            primaryDone = false;
            rm.simulateEvent(Track(parent, Vector3D(), Vector3D(0, 0, 1), kinetic), 1e30, 1000);
        }
        const double meanTau = sumTau / decays, meanL = sumL / decays;
        char label[80];
        std::snprintf(label, sizeof(label), "%s: proper time vs exp(-t/tau)", c.name);
        std::printf("  %-6s %8.2f | %10.4g %10.4g %+7.2f | %12.4g %12.4g %+7.2f\n", c.name, c.pGeV, meanTau, tau,
                    (meanTau - tau) / (tau / std::sqrt(decays)), meanL, meanLength,
                    (meanL - meanLength) / (meanLength / std::sqrt(decays)));
        printChi2(label, h, expectedCounts([](double t) { return std::exp(-t); }, 0.0, 5.0, nBins, decays));
    }
    std::printf("  pull = (MC - expected) / (sigma / sqrt(N)) for an exponential (sigma = mean).\n");
}

// ==============================================================================
// TEST 6  Decay law with energy loss and stopping
// ==============================================================================
// A pi+ / mu+ shot into a thick lead block slows down, stops, and decays at rest
// (or, rarely, in flight). Its proper time at the decay must STILL be exponential
// with mean tau: the clock runs at the rate 1/gamma while it moves, at rate 1 at
// rest, and the decay happens when it reaches the time drawn at the start.
void testDecayLawStopping(int n) {
    header("TEST 6  Decay law with energy loss and stopping (thick lead)");
    std::printf("  %-14s | %9s %9s | %10s %10s %7s\n", "primary", "in flight", "at rest", "<tau> [ns]", "PDG", "pull");
    struct Case { const char* name; double kineticMeV; };
    for (const Case& c : {Case{"pi+", 200.0}, Case{"mu+", 100.0}, Case{"pi-", 50.0}}) {
        const Particle* parent = find(c.name);
        const double tau = parent->getLifetime();
        auto setup = makeSetup(50.0, 50.0, 50.0);
        bool primaryDone = false;
        int inFlight = 0, atRest = 0;
        double sumTau = 0.0;
        const int nBins = 40;
        std::vector<double> h(nBins, 0.0);
        setup->runManager.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (primaryDone || step.processDefinedStep != "Decay") return;
            primaryDone = true;
            // At rest the zero-length at-rest step carries the decay
            if (step.stepLength > 0.0) ++inFlight; else ++atRest;
            sumTau += track.getProperTime();
            fill(h, track.getProperTime() / tau, 0.0, 5.0);
        });
        for (int i = 0; i < n; ++i) {
            primaryDone = false;
            setup->runManager.simulateEvent(Track(parent, Vector3D(0, 0, -49.0), Vector3D(0, 0, 1), c.kineticMeV * Units::MeV),
                                            1.0 * Units::m, 1000000);
        }
        const int decays = inFlight + atRest;
        const double meanTau = sumTau / decays;
        char label[80];
        std::snprintf(label, sizeof(label), "%s %.0f MeV: proper time vs exp(-t/tau)", c.name, c.kineticMeV);
        std::printf("  %-4s %5.0f MeV | %8.2f%% %8.2f%% | %10.4g %10.4g %+7.2f\n", c.name, c.kineticMeV,
                    100.0 * inFlight / decays, 100.0 * atRest / decays, meanTau, tau, (meanTau - tau) / (tau / std::sqrt(decays)));
        printChi2(label, h, expectedCounts([](double t) { return std::exp(-t); }, 0.0, 5.0, nBins, decays));
    }
    std::printf("  [INFO] A real pi- (and 97 %% of the mu- in lead) is captured by a nucleus at rest\n");
    std::printf("         instead of decaying: no strong interaction / muon capture here.\n");
}

// ==============================================================================
// TEST 7  Energy conservation in full events
// ==============================================================================
// Kinetic energy is NOT conserved: decays turn rest mass into kinetic energy,
// pair production does the opposite. With Q(decay) = m(parent) - sum m(daughters):
//   deposited + escaped + neutrinos = T0 + sum Q(decays) + 2 m_e (N_annih - N_conv)
// exactly. (Equivalently: total energy is conserved.) The neutrinos never
// interact, so their energy is exactly what leaves the detector unseen.
void testEnergyConservation(int nEvents) {
    header("TEST 7  Energy balance with decays (lead block 20 x 20 x 20 cm), fractions of T0 + Q");
    std::printf("%-15s %5s | %9s %9s %9s | %8s %7s %7s | %12s %7s\n", "primary", "evts", "deposited", "escaped",
                "neutrinos", "Q/(T0+Q)", "N_dec", "N_conv", "max |imbal.|", "tracks");
    struct Case { const char* name; double energy; const char* label; };
    for (const Case& c : {Case{"mu+", 0.1, "mu+ 100 MeV"}, Case{"mu-", 1.0, "mu- 1 GeV"}, Case{"pi+", 0.3, "pi+ 300 MeV"},
                          Case{"pi-", 1.0, "pi- 1 GeV"}, Case{"pi0", 1.0, "pi0 1 GeV"}, Case{"tau-", 10.0, "tau- 10 GeV"},
                          Case{"tau+", 1.0, "tau+ 1 GeV"}}) {
        const Particle* primaryDef = find(c.name);
        const double t0 = c.energy * Units::GeV;
        auto setup = makeSetup(10.0, 10.0, 10.0);
        double dep = 0.0, esc = 0.0, nu = 0.0, q = 0.0;
        int nConv = 0, nAnnih = 0, nDecay = 0;
        setup->runManager.setUserSteppingAction([&](const Track& track, const Step& step) {
            dep += step.totalEnergyDeposit;
            if (step.leftWorld) {
                if (isNeutrino(*track.getParticleDefinition())) nu += track.getKineticEnergy();
                else esc += track.getKineticEnergy();
            }
            if (step.processDefinedStep == "conv") ++nConv;
            if (step.processDefinedStep == "annihil") ++nAnnih;
            if (step.processDefinedStep == "Decay") {
                ++nDecay;
                double qStep = track.getMass();
                for (const Track& d : *step.secondaryInCurrentStep) qStep -= d.getMass();
                q += qStep;
            }
        });
        double sumDep = 0.0, sumEsc = 0.0, sumNu = 0.0, sumQ = 0.0, maxImbalance = 0.0;
        long sumDecay = 0, sumConv = 0, sumTracks = 0;
        for (int i = 0; i < nEvents; ++i) {
            dep = esc = nu = q = 0.0;
            nConv = nAnnih = nDecay = 0;
            setup->runManager.simulateEvent(Track(primaryDef, Vector3D(0, 0, -10.5), Vector3D(0, 0, 1), t0), 1.0 * Units::m, 1000000);
            const double available = t0 + q + 2.0 * m_e * (nAnnih - nConv);
            maxImbalance = std::max(maxImbalance, std::abs(dep + esc + nu - available) / (t0 + q));
            sumDep += dep / (t0 + q);
            sumEsc += esc / (t0 + q);
            sumNu += nu / (t0 + q);
            sumQ += q / (t0 + q);
            sumDecay += nDecay;
            sumConv += nConv;
            sumTracks += setup->runManager.getNumberOfTracks();
        }
        std::printf("%-15s %5d | %8.2f%% %8.2f%% %8.2f%% | %7.2f%% %7.2f %7.1f | %12.1e %7.0f\n", c.label, nEvents,
                    100.0 * sumDep / nEvents, 100.0 * sumEsc / nEvents, 100.0 * sumNu / nEvents, 100.0 * sumQ / nEvents,
                    static_cast<double>(sumDecay) / nEvents, static_cast<double>(sumConv) / nEvents, maxImbalance,
                    static_cast<double>(sumTracks) / nEvents);
    }
    std::printf("deposited / escaped / neutrinos: mean fractions of T0 + Q (Q = energy released by decays);\n");
    std::printf("they add up to 100 %% up to the 2 m_e terms. N_dec: decays per event. imbalance: rounding only.\n");
}

} // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::string(argv[1]) == "quick";
    RandomEngine::setSeed(12345);
    std::printf("[RandomEngine] Seed: %llu%s\n", static_cast<unsigned long long>(RandomEngine::getSeed()),
                quick ? "  (quick mode)" : "");

    const int f = quick ? 10 : 1;
    testMuonDecay(1000000 / f);
    testTauLeptonic(1000000 / f);
    testBranchingRatios(1000000 / f);
    testKinematics(2000 / f);
    testDecayLawVacuum(100000 / f);
    testDecayLawStopping(20000 / f);
    testEnergyConservation(200 / f);
    return 0;
}
