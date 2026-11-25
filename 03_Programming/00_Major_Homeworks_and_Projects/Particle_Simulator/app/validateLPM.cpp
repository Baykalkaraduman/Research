// ==============================================================================
// [STAGE 9] Validation of the Landau-Pomeranchuk-Migdal (LPM) and Ter-Mikaelian
//           dielectric suppression effects in lead.
//
//   TEST 1  Migdal functions phi(s), G(s), xi(s), solver residual & Pb scales
//   TEST 2  Bremsstrahlung suppression: 25 GeV SLAC E-146 regime & TeV-PeV dE/dx
//   TEST 3  Pair production suppression: XCOM match at 100 GeV, U-shape & PeV scaling
//   TEST 4  MC sampling fidelity: bremsstrahlung & pair quantile chi2/ndf + pulls
//   TEST 5  100 TeV gamma shower in Pb: first-conversion MFP & shower elongation
//
// Run:  ./ValidationLPMApp
//       ./ValidationLPMApp quick
//       ./ValidationLPMApp 42
// ==============================================================================

#include "Box.h"
#include "Bremsstrahlung.h"
#include "GammaConversion.h"
#include "Ionization.h"
#include "LPMFunctions.h"
#include "Material.h"
#include "Particle.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using PhysicalConstants::m_e;

namespace {

const Particle gammaDef(22, "gamma", 0.0, 0.0, 1.0, -1.0);
const Particle electronDef(11, "e-", m_e, -1.0, 0.5, -1.0);
const Particle positronDef(-11, "e+", m_e, 1.0, 0.5, -1.0);

constexpr double kMaxPull = 3.0;
constexpr double kMaxChi2PerNdf = 3.0;

std::shared_ptr<Material> makeLead() {
    const DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

void printHeader(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

struct QuantileRef {
    double mean = 0.0;
    double variance = 0.0;
    std::array<double, 9> edges{}; // 9 interior bin edges for 10 equal-probability bins
};

template <typename PdfFunc, typename ObsFunc>
QuantileRef buildQuantileRef(double uMin, double uMax, PdfFunc pdfInU, ObsFunc obsInU) {
    constexpr int kSteps = 2000;
    const double du = (uMax - uMin) / kSteps;
    std::vector<double> uGrid(static_cast<std::size_t>(kSteps + 1));
    std::vector<double> pGrid(static_cast<std::size_t>(kSteps + 1));
    std::vector<double> cdf(static_cast<std::size_t>(kSteps + 1), 0.0);

    for (int i = 0; i <= kSteps; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        uGrid[idx] = uMin + i * du;
        pGrid[idx] = std::max(0.0, pdfInU(uGrid[idx]));
        if (i > 0) {
            cdf[idx] = cdf[idx - 1] + 0.5 * (pGrid[idx - 1] + pGrid[idx]) * du;
        }
    }

    QuantileRef ref;
    const double total = cdf.back();
    if (total <= 0.0) {
        return ref;
    }

    double sumObs = 0.0;
    double sumObs2 = 0.0;
    for (int i = 1; i <= kSteps; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        const double dA = cdf[idx] - cdf[idx - 1];
        const double obs = obsInU(0.5 * (uGrid[idx - 1] + uGrid[idx]));
        sumObs += dA * obs;
        sumObs2 += dA * obs * obs;
    }
    ref.mean = sumObs / total;
    ref.variance = std::max(0.0, sumObs2 / total - ref.mean * ref.mean);

    for (int b = 1; b < 10; ++b) {
        const double target = (static_cast<double>(b) / 10.0) * total;
        const auto it = std::lower_bound(cdf.begin(), cdf.end(), target);
        if (it == cdf.begin()) {
            ref.edges[static_cast<std::size_t>(b - 1)] = uGrid.front();
        } else if (it == cdf.end()) {
            ref.edges[static_cast<std::size_t>(b - 1)] = uGrid.back();
        } else {
            const auto hi = static_cast<std::size_t>(it - cdf.begin());
            const auto lo = hi - 1;
            const double span = cdf[hi] - cdf[lo];
            const double frac = (span > 0.0) ? (target - cdf[lo]) / span : 0.0;
            ref.edges[static_cast<std::size_t>(b - 1)] = uGrid[lo] + frac * (uGrid[hi] - uGrid[lo]);
        }
    }
    return ref;
}

int findBin(const std::array<double, 9>& edges, double u) {
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (u < edges[i]) {
            return static_cast<int>(i);
        }
    }
    return static_cast<int>(edges.size());
}

// ==============================================================================
// TEST 1  Migdal functions & lead characteristic scales
// ==============================================================================
bool testMigdalFunctions() {
    printHeader("TEST 1  Migdal functions phi(s), G(s), xi(s) & Pb characteristic scales");
    auto lead = makeLead();
    bool ok = true;

    const double eLpmTeV = LPMFunctions::computeLPMEnergy(*lead) / Units::TeV;
    const double hbarWpEv = LPMFunctions::computePlasmaEnergy(*lead) / Units::eV;
    std::printf("Lead (X0 = %.4f cm): E_LPM = %.4f TeV (ref ~4.312 TeV), hbar*w_p = %.2f eV (ref ~61.1 eV)\n\n",
                lead->getRadLength() / Units::cm, eLpmTeV, hbarWpEv);

    if (std::abs(eLpmTeV - 4.3125) > 0.01 || std::abs(hbarWpEv - 61.1) > 0.5) {
        ok = false;
    }

    // Closed-form identities derived from cotangent / hyperbolic series
    const double phi05 = LPMFunctions::computePhi(0.5);
    const double g05 = LPMFunctions::computeG(0.5);
    const double phi10 = LPMFunctions::computePhi(1.0);
    const double g10 = LPMFunctions::computeG(1.0);

    constexpr double kRefPhi05 = 0.93156811;
    constexpr double kRefG05 = 0.87251456;
    constexpr double kRefPhi10 = 0.99062121;
    constexpr double kRefG10 = 0.98196274;

    std::printf("Exact points: phi(0.5)=%.6f (ref %.6f), G(0.5)=%.6f (ref %.6f)\n", phi05, kRefPhi05, g05, kRefG05);
    std::printf("              phi(1.0)=%.6f (ref %.6f), G(1.0)=%.6f (ref %.6f)\n", phi10, kRefPhi10, g10, kRefG10);

    if (std::abs(phi05 - kRefPhi05) > 1.0e-4 || std::abs(g05 - kRefG05) > 1.0e-4 ||
        std::abs(phi10 - kRefPhi10) > 1.0e-4 || std::abs(g10 - kRefG10) > 1.0e-4) {
        ok = false;
    }

    // Small-s limits: phi(s) -> 6s, G(s) -> 12*pi*s^2
    constexpr double sSmall = 1.0e-3;
    const double phiRatio = LPMFunctions::computePhi(sSmall) / (6.0 * sSmall);
    const double gRatio = LPMFunctions::computeG(sSmall) / (12.0 * M_PI * sSmall * sSmall);
    const double phiJump = std::abs(LPMFunctions::computePhi(1.999999) - LPMFunctions::computePhi(2.000001));
    const double gJump = std::abs(LPMFunctions::computeG(1.999999) - LPMFunctions::computeG(2.000001));

    std::printf("Small-s (1e-3): phi/(6s) = %.5f, G/(12*pi*s^2) = %.5f | Join jump at s=2: dPhi=%.2e, dG=%.2e\n",
                phiRatio, gRatio, phiJump, gJump);
    if (std::abs(phiRatio - 1.0) > 0.01 || std::abs(gRatio - 1.0) > 0.01 || phiJump > 1.0e-4 || gJump > 1.0e-4) {
        ok = false;
    }

    // Fixed-point solver check: s * sqrt(xi(s)) == s'
    const double s1 = LPMFunctions::computeS1(82);
    double maxResidual = 0.0;
    for (double sPrime : {1.0e-4, 1.0e-2, 0.1, 0.5, 0.99, 2.0}) {
        const double s = LPMFunctions::solveMigdalS(sPrime, s1);
        const double xi = LPMFunctions::computeXi(s, s1);
        maxResidual = std::max(maxResidual, std::abs(s * std::sqrt(xi) - sPrime));
    }
    std::printf("Solver max |s*sqrt(xi(s)) - s'| = %.2e\n", maxResidual);
    if (maxResidual > 1.0e-10) {
        ok = false;
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 2  Bremsstrahlung suppression: 25 GeV SLAC E-146 regime & TeV-PeV dE/dx
// ==============================================================================
bool testBremsstrahlungSuppression() {
    printHeader("TEST 2  Bremsstrahlung LPM & dielectric suppression in lead");
    auto lead = makeLead();
    bool ok = true;

    const double e25 = 25.0 * Units::GeV;
    std::printf("25 GeV e- in Pb (k_LPM ~ 145 MeV, k_p = gamma*hbar*w_p ~ 3.0 MeV):\n");
    std::printf("%10s %12s %12s\n", "k [MeV]", "R (LPM only)", "R (LPM + DE)");
    for (double kMeV : {0.5, 2.0, 5.0, 20.0, 100.0, 500.0, 5000.0}) {
        const double k = kMeV * Units::MeV;
        const double rLpm = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, k, *lead, false);
        const double rDe = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, k, *lead, true);
        std::printf("%10.1f %12.4f %12.4f\n", kMeV, rLpm, rDe);
    }

    const double rHigh = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, 5.0 * Units::GeV, *lead, true);
    const double rMidLpm = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, 5.0 * Units::MeV, *lead, false);
    const double rLowLpm = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, 0.5 * Units::MeV, *lead, false);
    const double rLowDe = LPMFunctions::computeBremsstrahlungFactor(e25 + m_e, 0.5 * Units::MeV, *lead, true);

    if (std::abs(rHigh - 1.0) > 0.05 || rMidLpm < 0.25 || rMidLpm > 0.65 || rLowDe > 0.05 || rLowDe >= 0.25 * rLowLpm) {
        ok = false;
    }

    std::printf("\nTotal radiative stopping power dE/dx (LPM+DE vs Bethe-Heitler):\n");
    std::printf("%12s %14s %14s %10s\n", "T [GeV]", "BH [GeV/cm]", "LPM [GeV/cm]", "LPM / BH");

    double prevRatio = 2.0;
    for (double eGeV : {10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0}) {
        const double t = eGeV * Units::GeV;
        const double dedxBh = Bremsstrahlung::computeRadiativeDEDX(t, *lead, false, false);
        const double dedxLpm = Bremsstrahlung::computeRadiativeDEDX(t, *lead, true, true);
        const double ratio = dedxLpm / dedxBh;
        std::printf("%12.0f %14.4e %14.4e %10.4f\n", eGeV, dedxBh / Units::GeV, dedxLpm / Units::GeV, ratio);

        if (eGeV <= 10.0 && std::abs(ratio - 1.0) > 0.005) {
            ok = false;
        }
        if (eGeV >= 1000.0 && ratio >= prevRatio) {
            ok = false;
        }
        if (eGeV >= 1000000.0 && ratio > 0.35) {
            ok = false;
        }
        prevRatio = ratio;
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 3  Pair production LPM suppression
// ==============================================================================
bool testPairProductionSuppression() {
    printHeader("TEST 3  Pair production LPM suppression in lead (100 GeV - 10 PeV)");
    auto lead = makeLead();
    bool ok = true;

    const double barn = 1.0e-24 * Units::cm * Units::cm;
    std::printf("%12s %12s %12s %10s %10s\n", "E [TeV]", "BH [barn]", "LPM [barn]", "LPM / BH", "R(x=0.5)");

    double sigma1PeV = 0.0;
    double sigma10PeV = 0.0;
    double prevRatio = 2.0;

    for (double eTeV : {0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0}) {
        const double e = eTeV * Units::TeV;
        const double sBh = GammaConversion::getCrossSectionPerAtom(e, *lead, false);
        const double sLpm = GammaConversion::getCrossSectionPerAtom(e, *lead, true);
        const double ratio = sLpm / sBh;
        const double rMid = LPMFunctions::computePairProductionFactor(e, 0.5, *lead);

        std::printf("%12.1f %12.4f %12.4f %10.4f %10.4f\n", eTeV, sBh / barn, sLpm / barn, ratio, rMid);

        if (eTeV <= 0.1 && std::abs(ratio - 1.0) > 0.002) {
            ok = false;
        }
        if (eTeV >= 10.0 && ratio >= prevRatio) {
            ok = false;
        }
        if (std::abs(eTeV - 1000.0) < 1.0) sigma1PeV = sLpm;
        if (std::abs(eTeV - 10000.0) < 1.0) sigma10PeV = sLpm;
        prevRatio = ratio;
    }

    // U-shaped differential spectrum at 100 TeV: central x=0.5 strongly dipped below wing x=0.05
    const double e100 = 100.0 * Units::TeV;
    const double dMid = GammaConversion::computeDifferentialCrossSection(e100, 0.5, *lead, true);
    const double dWing = GammaConversion::computeDifferentialCrossSection(e100, 0.05, *lead, true);
    const double peVScaleRatio = sigma10PeV / sigma1PeV; // ~ 1/sqrt(10) = 0.316 (plus mild log xi)

    std::printf("\n100 TeV pair dip: dsigma(0.5)/dsigma(0.05) = %.4f (must be < 0.70)\n", dMid / dWing);
    std::printf("PeV scaling: sigma(10 PeV) / sigma(1 PeV) = %.4f (1/sqrt(10) = 0.316)\n", peVScaleRatio);

    if (dMid / dWing >= 0.70 || peVScaleRatio < 0.25 || peVScaleRatio > 0.42) {
        ok = false;
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 4  MC sampling fidelity (Bremsstrahlung & Pair Production with LPM)
// ==============================================================================
bool testSamplingFidelity(int nSamples) {
    printHeader("TEST 4  MC sampling fidelity with LPM (10 quantile bins + mean pull)");
    auto lead = makeLead();
    bool ok = true;

    // 4a. Bremsstrahlung k sampling in u = ln(k)
    std::printf("Bremsstrahlung samplePhotonEnergy (LPM+DE, %d samples/energy):\n", nSamples);
    std::printf("%10s %10s %12s %12s %8s %10s\n", "T [GeV]", "kMin[MeV]", "<ln k> ref", "<ln k> MC", "pull", "chi2/ndf");

    struct BremCase {
        double tGeV;
        double kMinMeV;
    };
    for (const auto& c : std::array<BremCase, 2>{{{25.0, 1.0}, {100000.0, 1000.0}}}) {
        const double t = c.tGeV * Units::GeV;
        const double kMin = c.kMinMeV * Units::MeV;
        const double uMin = std::log(kMin);
        const double uMax = std::log(t);

        const auto ref = buildQuantileRef(
            uMin, uMax,
            [&](double u) {
                const double k = std::exp(u);
                return k * Bremsstrahlung::computeDifferentialCrossSection(t, k, *lead, true, true);
            },
            [](double u) {return u;});

        std::array<long, 10> counts{};
        double sumU = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            const double k = Bremsstrahlung::samplePhotonEnergy(t, kMin, *lead, true, true);
            const double u = std::log(k);
            sumU += u;
            ++counts[static_cast<std::size_t>(findBin(ref.edges, u))];
        }

        const double meanMc = sumU / nSamples;
        const double pull = (meanMc - ref.mean) / std::sqrt(ref.variance / nSamples);
        const double expBin = static_cast<double>(nSamples) / 10.0;
        double chi2 = 0.0;
        for (long count : counts) {
            const double diff = static_cast<double>(count) - expBin;
            chi2 += (diff * diff) / expBin;
        }
        const double chi2PerNdf = chi2 / 9.0;

        if (std::abs(pull) >= kMaxPull || chi2PerNdf >= kMaxChi2PerNdf) {
            ok = false;
        }
        std::printf("%10.0f %10.1f %12.5f %12.5f %+8.2f %10.2f\n",
                    c.tGeV, c.kMinMeV, ref.mean, meanMc, pull, chi2PerNdf);
    }

    // 4b. Pair production eps sampling in u = ln(eps) on [eps0, 0.5]
    std::printf("\nPair production sampleEnergyFraction (LPM, %d samples/energy):\n", nSamples);
    std::printf("%10s %10s %10s %10s %8s %10s\n", "E [TeV]", "<e> BH", "<e> LPM", "<e> MC", "pull", "chi2/ndf");

    for (double eTeV : {10.0, 100.0}) {
        const double e = eTeV * Units::TeV;
        const double epsMin = (m_e / e) * (1.0 + 1.0e-9);
        const double uMin = std::log(epsMin);
        const double uMax = std::log(0.5);

        const auto bhRef = buildQuantileRef(
            uMin, uMax,
            [&](double u) {
                const double eps = std::exp(u);
                return eps * GammaConversion::computeDifferentialCrossSection(e, eps, *lead, false);
            },
            [](double u) {return std::exp(u);});

        const auto lpmRef = buildQuantileRef(
            uMin, uMax,
            [&](double u) {
                const double eps = std::exp(u);
                return eps * GammaConversion::computeDifferentialCrossSection(e, eps, *lead, true);
            },
            [](double u) {return std::exp(u);});

        std::array<long, 10> counts{};
        double sumEps = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            const double eps = GammaConversion::sampleEnergyFraction(e, *lead, true);
            sumEps += eps;
            ++counts[static_cast<std::size_t>(findBin(lpmRef.edges, std::log(eps)))];
        }

        const double meanMc = sumEps / nSamples;
        const double pull = (meanMc - lpmRef.mean) / std::sqrt(lpmRef.variance / nSamples);
        const double expBin = static_cast<double>(nSamples) / 10.0;
        double chi2 = 0.0;
        for (long count : counts) {
            const double diff = static_cast<double>(count) - expBin;
            chi2 += (diff * diff) / expBin;
        }
        const double chi2PerNdf = chi2 / 9.0;

        if (std::abs(pull) >= kMaxPull || chi2PerNdf >= kMaxChi2PerNdf || lpmRef.mean >= bhRef.mean) {
            ok = false;
        }
        std::printf("%10.0f %10.5f %10.5f %10.5f %+8.2f %10.2f\n",
                    eTeV, bhRef.mean, lpmRef.mean, meanMc, pull, chi2PerNdf);
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 5  100 TeV gamma shower elongation in lead (LPM vs Bethe-Heitler)
// ==============================================================================
bool testShowerElongation(int nShowers) {
    printHeader("TEST 5  100 TeV gamma shower in lead: conversion MFP & shower elongation");
    auto lead = makeLead();
    const double x0 = lead->getRadLength();
    const double e0 = 100.0 * Units::TeV;
    const double showerCut = 200.0 * Units::GeV; // keeps core shower fast while deep in LPM regime

    struct ShowerResult {
        double meanConvX0 = 0.0;
        double refConvX0 = 0.0;
        double convPull = 0.0;
        double meanDepthX0 = 0.0;
        double errDepthX0 = 0.0;
        double maxRelEnergyErr = 0.0;
    };

    const auto runBatch = [&](bool enableLPM) -> ShowerResult {
        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        const double halfZ = 50.0 * x0; // 100 X0 thick block: z in [-50 X0, +50 X0]
        rm.setWorldVolume(std::make_unique<Box>("World", 50.0 * Units::cm, 50.0 * Units::cm, halfZ + 10.0 * Units::cm, vacuum));
        rm.addVolume(std::make_unique<Box>("Lead", 40.0 * Units::cm, 40.0 * Units::cm, halfZ, lead));

        rm.addPhysicsProcess(std::make_unique<Ionization>());
        auto brem = std::make_unique<Bremsstrahlung>(&gammaDef);
        brem->setProductionCut(showerCut);
        brem->setLPM(enableLPM);
        rm.addPhysicsProcess(std::move(brem));

        auto conv = std::make_unique<GammaConversion>(&electronDef, &positronDef);
        conv->setLPM(enableLPM);
        rm.addPhysicsProcess(std::move(conv));

        rm.setEnergyCut(showerCut);
        rm.setVerboseLevel(0);

        const double muConv = lead->getAtomDensity() * GammaConversion::getCrossSectionPerAtom(e0, *lead, enableLPM);
        const double refLambdaX0 = (1.0 / muConv) / x0;

        bool converted = false;
        double convDepthSum = 0.0;
        long nConverted = 0;

        double depEnergyEvent = 0.0;
        double depZWeightEvent = 0.0;
        double escapedEvent = 0.0;

        double sumDepth = 0.0;
        double sumDepth2 = 0.0;
        double maxRelErr = 0.0;

        rm.setUserSteppingAction([&](const Track& track, const Step& step) {
            const double zDepthX0 = (step.postStepPosition.z + halfZ) / x0;
            if (!converted && step.processDefinedStep == "conv") {
                converted = true;
                convDepthSum += zDepthX0;
                ++nConverted;
            }
            if (step.totalEnergyDeposit > 0.0) {
                const double zMidX0 = (0.5 * (step.preStepPosition.z + step.postStepPosition.z) + halfZ) / x0;
                depEnergyEvent += step.totalEnergyDeposit;
                depZWeightEvent += step.totalEnergyDeposit * zMidX0;
            }
            if (step.leftWorld) {
                escapedEvent += track.getKineticEnergy();
            }
        });

        for (int i = 0; i < nShowers; ++i) {
            converted = false;
            depEnergyEvent = 0.0;
            depZWeightEvent = 0.0;
            escapedEvent = 0.0;

            Track gamma(&gammaDef, Vector3D(0.0, 0.0, -halfZ + 1.0e-5), Vector3D(0.0, 0.0, 1.0), e0);
            rm.simulateEvent(gamma, 2.0 * Units::cm, 10000);

            if (depEnergyEvent > 0.0) {
                const double depthX0 = depZWeightEvent / depEnergyEvent;
                sumDepth += depthX0;
                sumDepth2 += depthX0 * depthX0;
            }
            maxRelErr = std::max(maxRelErr, std::abs((depEnergyEvent + escapedEvent) - e0) / e0);
        }

        ShowerResult r;
        r.refConvX0 = refLambdaX0;
        r.meanConvX0 = (nConverted > 0) ? convDepthSum / static_cast<double>(nConverted) : 0.0;
        r.convPull = (nConverted > 0)
                   ? (r.meanConvX0 - refLambdaX0) / (refLambdaX0 / std::sqrt(static_cast<double>(nConverted)))
                   : 0.0;
        r.meanDepthX0 = sumDepth / nShowers;
        const double varDepth = std::max(0.0, sumDepth2 / nShowers - r.meanDepthX0 * r.meanDepthX0);
        r.errDepthX0 = std::sqrt(varDepth / nShowers);
        r.maxRelEnergyErr = maxRelErr;
        return r;
    };

    const ShowerResult bh = runBatch(false);
    const ShowerResult lpm = runBatch(true);

    std::printf("(%d showers/mode, 100 X0 Pb block, core cut = 200 GeV)\n", nShowers);
    std::printf("%-8s | %10s %10s %7s | %14s | %10s\n",
                "Mode", "<z_conv>", "lambda", "pull", "<z_dep> [X0]", "max|dE|/E");
    std::printf("%-8s | %8.2f X0 %8.2f X0 %+7.2f | %6.2f +/- %.2f | %10.2e\n",
                "BH", bh.meanConvX0, bh.refConvX0, bh.convPull, bh.meanDepthX0, bh.errDepthX0, bh.maxRelEnergyErr);
    std::printf("%-8s | %8.2f X0 %8.2f X0 %+7.2f | %6.2f +/- %.2f | %10.2e\n",
                "LPM+DE", lpm.meanConvX0, lpm.refConvX0, lpm.convPull, lpm.meanDepthX0, lpm.errDepthX0, lpm.maxRelEnergyErr);

    const double shiftX0 = lpm.meanDepthX0 - bh.meanDepthX0;
    const double shiftSigma = shiftX0 / std::sqrt(bh.errDepthX0 * bh.errDepthX0 + lpm.errDepthX0 * lpm.errDepthX0);
    std::printf("LPM longitudinal shift: %+.2f X0 (%+.1f sigma)\n", shiftX0, shiftSigma);

    const bool ok = (std::abs(bh.convPull) < kMaxPull) &&
                    (std::abs(lpm.convPull) < kMaxPull) &&
                    (shiftX0 > 3.0) && (shiftSigma > 5.0) &&
                    (bh.maxRelEnergyErr < 1.0e-3) && (lpm.maxRelEnergyErr < 1.0e-3);

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

} // namespace

int main(int argc, char** argv) {
    bool quick = false;
    unsigned long long seed = 12345;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "quick") quick = true;
        else seed = std::stoull(arg);
    }
    RandomEngine::setSeed(seed);

    std::printf("[RandomEngine] Seed: %llu%s\n", static_cast<unsigned long long>(RandomEngine::getSeed()),
                quick ? "   (quick mode)" : "");

    bool allPassed = true;
    allPassed &= testMigdalFunctions();
    allPassed &= testBremsstrahlungSuppression();
    allPassed &= testPairProductionSuppression();
    allPassed &= testSamplingFidelity(quick ? 10000 : 100000);
    allPassed &= testShowerElongation(quick ? 100 : 400);

    std::printf("\n==============================================================================\n");
    std::printf("  OVERALL: %s\n", allPassed ? "PASS" : "FAIL");
    std::printf("==============================================================================\n");
    return allPassed ? 0 : 1;
}