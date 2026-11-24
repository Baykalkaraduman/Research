// ==============================================================================
// Validation of Rayleigh scattering and bound Compton scattering in lead.
//
// Data: ENDF/B-VIII.0 (EPICS2017) in data/photoat-082_Pb_000.endf
//       (MF=23 MT=502/504 cross sections, MF=27 MT=502/504 F(x) and S(x))
//       and NIST XCOM in data/xcom_Z82.txt.
//
//   TEST 1  F(x) & S(x) limits and numerical integral of dsigma vs MF=23 MT=502/504
//   TEST 2  ENDF MF=23 MT=502 (coherent) & MT=504 (incoherent) vs NIST XCOM
//   TEST 3  Rayleigh angular sampling vs (1 + cos^2 theta) F^2(x) (quantile chi2 + pull)
//   TEST 4  Bound Compton sampling vs dsigma_KN/deps * S(x)/Z (forward suppression + chi2)
//   TEST 5  Event-by-event energy conservation & thin-slab transmission (phot+compt+Rayl)
//
// Run:  ./ValidationScatteringApp
//       ./ValidationScatteringApp quick
//       ./ValidationScatteringApp 42
// ==============================================================================

#include "Box.h"
#include "ComptonScattering.h"
#include "Ionization.h"
#include "Material.h"
#include "Particle.h"
#include "PhotoElectricEffect.h"
#include "PhotoElectricShellData.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RayleighScattering.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

namespace {

const Particle gammaDef(22, "gamma", 0.0, 0.0, 1.0, -1.0);
const Particle electronDef(11, "e-", m_e, -1.0, 0.5, -1.0);

constexpr int kZ = 82;
constexpr double kHcEvAngstrom = 12398.419843320026; // h*c [eV * Angstrom]
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

struct DistributionSummary {
    double integral = 0.0;      // [cm^2]
    double meanObservable = 0.0;
    double varObservable = 0.0;
    double forwardFraction = 0.0;
    std::array<double, 9> quantileX{}; // 9 interior bin edges in x [1/Angstrom] for 10 equal-prob bins
};

// Integrates a differential cross section on x in [0, xMax] using stretched variable
// x(t) = x0 * (exp(t) - 1) so both small-angle peaks and large-angle tails are resolved.
template <typename PdfFunc, typename ObsFunc>
DistributionSummary integrateInX(double xMax, PdfFunc pdfInX, ObsFunc observableAtX) {
    constexpr int kSteps = 2000;
    constexpr double kX0 = 0.05; // [1/Angstrom]
    const double tMax = std::log1p(xMax / kX0);
    const double dt = tMax / kSteps;

    std::vector<double> xGrid(static_cast<std::size_t>(kSteps + 1));
    std::vector<double> wGrid(static_cast<std::size_t>(kSteps + 1)); // pdf(x) * dx/dt
    std::vector<double> cdf(static_cast<std::size_t>(kSteps + 1), 0.0);

    for (int i = 0; i <= kSteps; ++i) {
        const double t = i * dt;
        const double expT = std::exp(t);
        const double x = std::min(xMax, kX0 * (expT - 1.0));
        const double dxdt = kX0 * expT;
        xGrid[static_cast<std::size_t>(i)] = x;
        wGrid[static_cast<std::size_t>(i)] = pdfInX(x) * dxdt;
        if (i > 0) {
            const auto idx = static_cast<std::size_t>(i);
            cdf[idx] = cdf[idx - 1] + 0.5 * (wGrid[idx - 1] + wGrid[idx]) * dt;
        }
    }

    DistributionSummary out;
    out.integral = cdf.back();
    if (out.integral <= 0.0) {
        return out;
    }

    double sumObs = 0.0;
    double sumObs2 = 0.0;
    double sumForward = 0.0;
    for (int i = 1; i <= kSteps; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        const double dA = cdf[idx] - cdf[idx - 1];
        const double xMid = 0.5 * (xGrid[idx - 1] + xGrid[idx]);
        const double obs = observableAtX(xMid);
        sumObs += dA * obs;
        sumObs2 += dA * obs * obs;
        const double cosTheta = 1.0 - 2.0 * (xMid / xMax) * (xMid / xMax);
        if (cosTheta > 0.5) {
            sumForward += dA;
        }
    }

    out.meanObservable = sumObs / out.integral;
    out.varObservable = std::max(0.0, sumObs2 / out.integral - out.meanObservable * out.meanObservable);
    out.forwardFraction = sumForward / out.integral;

    constexpr int kBins = 10;
    for (int b = 1; b < kBins; ++b) {
        const double target = (static_cast<double>(b) / kBins) * out.integral;
        const auto it = std::lower_bound(cdf.begin(), cdf.end(), target);
        if (it == cdf.begin()) {
            out.quantileX[static_cast<std::size_t>(b - 1)] = xGrid.front();
        } else if (it == cdf.end()) {
            out.quantileX[static_cast<std::size_t>(b - 1)] = xGrid.back();
        } else {
            const std::size_t hi = static_cast<std::size_t>(it - cdf.begin());
            const std::size_t lo = hi - 1;
            const double span = cdf[hi] - cdf[lo];
            const double frac = (span > 0.0) ? (target - cdf[lo]) / span : 0.0;
            out.quantileX[static_cast<std::size_t>(b - 1)] = xGrid[lo] + frac * (xGrid[hi] - xGrid[lo]);
        }
    }
    return out;
}

int findQuantileBin(const std::array<double, 9>& edges, double x) {
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (x < edges[i]) {
            return static_cast<int>(i);
        }
    }
    return static_cast<int>(edges.size());
}

// ==============================================================================
// TEST 1  F(x), S(x), f'(E), f''(E) and differential cross-section integrals
// ==============================================================================
bool testScatteringFactorsAndIntegrals() {
    printHeader("TEST 1  F(x), S(x), f'/f'' and integral of dsigma vs ENDF MF=23");
    const PhotoElectricShellData& data = PhotoElectricShellData::get(kZ);
    const double barn = 1.0e-24 * Units::cm * Units::cm;
    const double hc = kHcEvAngstrom * Units::eV;

    bool ok = true;

    const double f0 = data.getFormFactor(0.0);
    const double fHigh = data.getFormFactor(1.0e3);
    const double s0 = data.getScatteringFunction(0.0);
    const double sHigh = data.getScatteringFunction(1.0e3);

    std::printf("F(0) = %.4f (Z = %d), F(1e3 A^-1) = %.3e\n", f0, kZ, fHigh);
    std::printf("S(0) = %.4f,          S(1e3 A^-1) = %.4f (Z = %d)\n\n", s0, sHigh, kZ);

    if (std::abs(f0 - kZ) > 1.0e-2 || fHigh > 0.1 || std::abs(s0) > 1.0e-6 || std::abs(sHigh / kZ - 1.0) > 1.0e-3) {
        ok = false;
    }

    std::printf("%8s | %7s %6s | %10s %10s %8s | %10s %10s %8s\n",
                "E [keV]", "f'", "f''", "Coh MT502", "Int(anom)", "diff", "Inc MT504", "Int(KN*S)", "diff");
    for (double eKeV : {10.0, 14.0, 50.0, 90.0, 200.0, 1000.0, 10000.0}) {
        const double e = eKeV * Units::keV;
        const double xMax = e / hc;
        const double k = e / m_e;
        const double fp = data.getRealAnomalousFactor(e);
        const double fpp = data.getImaginaryAnomalousFactor(e);
        const double fpp2 = fpp * fpp;

        const auto raylInt = integrateInX(
            xMax,
            [&](double x) {
                const double s = x / xMax;
                const double cosT = 1.0 - 2.0 * s * s;
                const double fReal = data.getFormFactor(x) + fp;
                const double modF2 = fReal * fReal + fpp2;
                return (4.0 * M_PI * r_e * r_e / (xMax * xMax)) * x * (1.0 + cosT * cosT) * modF2;
            },
            [&](double x) {
                const double s = x / xMax;
                return 1.0 - 2.0 * s * s;
            });

        const auto comptInt = integrateInX(
            xMax,
            [&](double x) {
                const double s = x / xMax;
                const double w = 2.0 * s * s;
                const double eps = 1.0 / (1.0 + k * w);
                const double dEpsDx = k * eps * eps * (4.0 * x / (xMax * xMax));
                return ComptonScattering::computeDifferentialCrossSection(e, eps) * dEpsDx * data.getScatteringFunction(x);
            },
            [&](double x) {
                const double s = x / xMax;
                return 1.0 / (1.0 + 2.0 * k * s * s);
            });

        const double cohEndf = data.getCoherentCrossSection(e);
        const double cohCalc = raylInt.integral / barn;
        const double cohDiff = cohCalc / cohEndf - 1.0;

        const double incEndf = data.getIncoherentCrossSection(e);
        const double incCalc = comptInt.integral / barn;
        const double incDiff = incCalc / incEndf - 1.0;

        if (std::abs(cohDiff) >= 0.02 || std::abs(incDiff) >= 0.02) {
            ok = false;
        }

        std::printf("%8.0f | %+7.3f %6.3f | %10.4f %10.4f %+7.2f%% | %10.4f %10.4f %+7.2f%%\n",
                    eKeV, fp, fpp, cohEndf, cohCalc, 100.0 * cohDiff, incEndf, incCalc, 100.0 * incDiff);
    }

    std::printf("Int(anom) integrates (1 + cos^2 theta) * [(F(x) + f'(E))^2 + f''(E)^2].\n");
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 2  ENDF MF=23 MT=502 & MT=504 vs NIST XCOM (xcom_Z82.txt)
// ==============================================================================
bool testCrossSectionsVsXcom() {
    printHeader("TEST 2  ENDF MF=23 coherent & incoherent vs NIST XCOM (Lead)");
    auto lead = makeLead();
    const RayleighScattering rayl;
    const ComptonScattering compt(&electronDef);

    std::ifstream in(std::string(TOYSIM_DATA_DIR) + "/xcom_Z82.txt");
    if (!in) {
        throw std::runtime_error("cannot open xcom_Z82.txt");
    }

    struct XcomRow {
        double energy; // [GeV]
        double coh;    // [cm^2/g]
        double inc;    // [cm^2/g]
    };

    std::vector<XcomRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double eMeV = 0.0, coh = 0.0, inc = 0.0, pe = 0.0, pp = 0.0;
        if (ss >> eMeV >> coh >> inc >> pe >> pp) {
            const double e = eMeV * Units::MeV;
            if (!rows.empty() && e == rows.back().energy) {
                continue; // skip duplicate absorption-edge row
            }
            rows.push_back({e, coh, inc});
        }
    }

    const double atomsPerGram = lead->getAtomDensity() / lead->getDensity();
    double maxIncMid = 0.0;
    double maxIncMidKeV = 0.0;
    double maxIncHigh = 0.0;
    double maxCohHigh = 0.0;

    for (const auto& row : rows) {
        const double eKeV = row.energy / Units::keV;
        const double incEndf = atomsPerGram * compt.computeCrossSectionPerAtom(row.energy, *lead);
        const double cohEndf = rayl.computeMassAttenuation(row.energy, *lead);

        if (row.inc > 0.0) {
            const double rel = std::abs(incEndf / row.inc - 1.0);
            if (eKeV >= 2.5 && eKeV <= 1.0e5) {
                if (rel > maxIncMid) {
                    maxIncMid = rel;
                    maxIncMidKeV = eKeV;
                }
            } else if (eKeV >= 6.0e5) {
                maxIncHigh = std::max(maxIncHigh, rel);
            }
        }
        if (row.coh > 0.0 && eKeV >= 100.0 && eKeV <= 1.0e5) {
            maxCohHigh = std::max(maxCohHigh, std::abs(cohEndf / row.coh - 1.0));
        }
    }

    std::printf("%11s | %10s %10s %8s | %10s %10s %8s\n",
                "E [keV]", "XCOM coh", "ENDF coh", "diff", "XCOM inc", "ENDF inc", "diff");
    for (double targetKeV : {3.0, 10.0, 50.0, 100.0, 300.0, 1000.0, 10000.0, 100000.0, 1000000.0}) {
        const double targetE = targetKeV * Units::keV;
        const auto it = std::min_element(rows.begin(), rows.end(), [&](const XcomRow& a, const XcomRow& b) {
            return std::abs(a.energy - targetE) < std::abs(b.energy - targetE);
        });
        const double cohEndf = rayl.computeMassAttenuation(it->energy, *lead);
        const double incEndf = atomsPerGram * compt.computeCrossSectionPerAtom(it->energy, *lead);
        const double cohDiff = (it->coh > 0.0) ? 100.0 * (cohEndf / it->coh - 1.0) : 0.0;
        const double incDiff = (it->inc > 0.0) ? 100.0 * (incEndf / it->inc - 1.0) : 0.0;
        std::printf("%11.1f | %10.4e %10.4e %+7.2f%% | %10.4e %10.4e %+7.2f%%\n",
                    it->energy / Units::keV, it->coh, cohEndf, cohDiff, it->inc, incEndf, incDiff);
    }

    const bool ok = (maxIncMid < 0.015);
    std::printf("\nIncoherent (MT=504) max |ENDF/XCOM - 1| in [2.5 keV, 100 MeV]: %.2f%% (at %.1f keV, limit 1.50%%)\n",
                100.0 * maxIncMid, maxIncMidKeV);
    std::printf("Incoherent (MT=504) max |ENDF/XCOM - 1| above 600 MeV (info): %.2f%%\n", 100.0 * maxIncHigh);
    std::printf("Coherent   (MT=502) max |ENDF/XCOM - 1| in [100 keV, 100 MeV] (info): %.2f%%\n", 100.0 * maxCohHigh);
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 3  Rayleigh angular sampling vs (1 + cos^2 theta) * F^2(x)
// ==============================================================================
bool testRayleighSampling(int nSamples) {
    printHeader("TEST 3  Rayleigh angular sampling vs (1 + cos^2 theta) * F^2(x)");
    const PhotoElectricShellData& data = PhotoElectricShellData::get(kZ);
    const RayleighScattering rayl;
    const double hc = kHcEvAngstrom * Units::eV;

    bool ok = true;
    std::printf("(%d samples per energy, 10 equal-probability quantile bins in x)\n", nSamples);
    std::printf("%10s %12s %12s %8s %10s\n", "E [keV]", "<cosT> ref", "<cosT> MC", "pull", "chi2/ndf");

    for (double eKeV : {10.0, 50.0, 200.0, 1000.0}) {
        const double e = eKeV * Units::keV;
        const double xMax = e / hc;

        const auto ref = integrateInX(
            xMax,
            [&](double x) {
                const double s = x / xMax;
                const double cosT = 1.0 - 2.0 * s * s;
                const double f = data.getFormFactor(x);
                return (4.0 * M_PI * r_e * r_e / (xMax * xMax)) * x * (1.0 + cosT * cosT) * f * f;
            },
            [&](double x) {
                const double s = x / xMax;
                return 1.0 - 2.0 * s * s;
            });

        std::array<long, 10> counts{};
        double sumCos = 0.0;
        for (int i = 0; i < nSamples; ++i) {
            const double cosT = rayl.sampleCosTheta(e, kZ);
            sumCos += cosT;
            const double x = xMax * std::sqrt(std::max(0.0, 0.5 * (1.0 - cosT)));
            ++counts[static_cast<std::size_t>(findQuantileBin(ref.quantileX, x))];
        }

        const double meanMc = sumCos / nSamples;
        const double sigmaMean = std::sqrt(ref.varObservable / nSamples);
        const double pull = (sigmaMean > 0.0) ? (meanMc - ref.meanObservable) / sigmaMean : 0.0;

        const double expectedPerBin = static_cast<double>(nSamples) / 10.0;
        double chi2 = 0.0;
        for (long c : counts) {
            const double diff = static_cast<double>(c) - expectedPerBin;
            chi2 += (diff * diff) / expectedPerBin;
        }
        const double chi2PerNdf = chi2 / 9.0;

        if (std::abs(pull) >= kMaxPull || chi2PerNdf >= kMaxChi2PerNdf) {
            ok = false;
        }
        std::printf("%10.0f %12.6f %12.6f %+8.2f %10.2f\n", eKeV, ref.meanObservable, meanMc, pull, chi2PerNdf);
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 4  Bound Compton sampling vs dsigma_KN/deps * S(x)/Z
// ==============================================================================
bool testBoundComptonSampling(int nSamples) {
    printHeader("TEST 4  Bound Compton sampling vs dsigma_KN/deps * S(x)/Z");
    const PhotoElectricShellData& data = PhotoElectricShellData::get(kZ);
    const double hc = kHcEvAngstrom * Units::eV;

    bool ok = true;
    std::printf("(%d samples per energy, 10 quantile bins; fwd = fraction with cos(theta) > 0.5)\n", nSamples);
    std::printf("%8s | %9s %9s %9s %7s | %8s %8s %8s %7s | %8s\n",
                "E[keV]", "<e> free", "<e> bnd", "<e> MC", "pull", "fwd free", "fwd bnd", "fwd MC", "pull", "chi2/ndf");

    for (double eKeV : {10.0, 50.0, 200.0, 1000.0}) {
        const double e = eKeV * Units::keV;
        const double xMax = e / hc;
        const double k = e / m_e;

        const auto freeRef = integrateInX(
            xMax,
            [&](double x) {
                const double s = x / xMax;
                const double eps = 1.0 / (1.0 + 2.0 * k * s * s);
                const double dEpsDx = k * eps * eps * (4.0 * x / (xMax * xMax));
                return ComptonScattering::computeDifferentialCrossSection(e, eps) * dEpsDx;
            },
            [&](double x) {
                const double s = x / xMax;
                return 1.0 / (1.0 + 2.0 * k * s * s);
            });

        const auto boundRef = integrateInX(
            xMax,
            [&](double x) {
                const double s = x / xMax;
                const double eps = 1.0 / (1.0 + 2.0 * k * s * s);
                const double dEpsDx = k * eps * eps * (4.0 * x / (xMax * xMax));
                return ComptonScattering::computeDifferentialCrossSection(e, eps) * dEpsDx * data.getScatteringFunction(x);
            },
            [&](double x) {
                const double s = x / xMax;
                return 1.0 / (1.0 + 2.0 * k * s * s);
            });

        std::array<long, 10> counts{};
        double sumEps = 0.0;
        long forwardCount = 0;
        for (int i = 0; i < nSamples; ++i) {
            const double eps = ComptonScattering::sampleEnergyFraction(e, kZ);
            sumEps += eps;
            const double oneMinusCos = std::clamp((1.0 - eps) / (eps * k), 0.0, 2.0);
            const double cosT = 1.0 - oneMinusCos;
            if (cosT > 0.5) {
                ++forwardCount;
            }
            const double x = xMax * std::sqrt(0.5 * oneMinusCos);
            ++counts[static_cast<std::size_t>(findQuantileBin(boundRef.quantileX, x))];
        }

        const double meanMc = sumEps / nSamples;
        const double sigmaEps = std::sqrt(boundRef.varObservable / nSamples);
        const double epsPull = (sigmaEps > 0.0) ? (meanMc - boundRef.meanObservable) / sigmaEps : 0.0;

        const double fwdMc = static_cast<double>(forwardCount) / nSamples;
        const double pFwd = boundRef.forwardFraction;
        const double sigmaFwd = std::sqrt(std::max(1.0e-12, pFwd * (1.0 - pFwd)) / nSamples);
        const double fwdPull = (fwdMc - pFwd) / sigmaFwd;

        const double expectedPerBin = static_cast<double>(nSamples) / 10.0;
        double chi2 = 0.0;
        for (long c : counts) {
            const double diff = static_cast<double>(c) - expectedPerBin;
            chi2 += (diff * diff) / expectedPerBin;
        }
        const double chi2PerNdf = chi2 / 9.0;

        if (std::abs(epsPull) >= kMaxPull || std::abs(fwdPull) >= kMaxPull || chi2PerNdf >= kMaxChi2PerNdf) {
            ok = false;
        }

        std::printf("%8.0f | %9.5f %9.5f %9.5f %+7.2f | %8.4f %8.4f %8.4f %+7.2f | %8.2f\n",
                    eKeV, freeRef.meanObservable, boundRef.meanObservable, meanMc, epsPull,
                    freeRef.forwardFraction, pFwd, fwdMc, fwdPull, chi2PerNdf);
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 5  Energy conservation & thin-slab transmission (phot + compt + Rayl)
// ==============================================================================
bool testConservationAndTransmission(int nEvents, int nPhotons) {
    printHeader("TEST 5  Energy conservation & thin-slab transmission (phot + compt + Rayl)");
    auto lead = makeLead();
    RayleighScattering rayl;
    ComptonScattering compt(&electronDef);
    PhotoElectricEffect phot(&electronDef);
    rayl.setVerboseLevel(0);
    compt.setVerboseLevel(0);
    phot.setVerboseLevel(0);

    bool ok = true;
    double maxRaylErr = 0.0;
    double maxComptErr = 0.0;
    double maxDirErr = 0.0;
    long raylSecondaries = 0;

    std::vector<Track> secondaries;
    secondaries.reserve(4);

    for (double eKeV : {10.0, 50.0, 200.0, 1000.0}) {
        const double e = eKeV * Units::keV;
        for (int n = 0; n < nEvents; ++n) {
            // Rayleigh postStepDoIt
            {
                Track gamma(&gammaDef, Vector3D(), Vector3D(0.0, 0.0, 1.0), e);
                Step step;
                secondaries.clear();
                rayl.postStepDoIt(gamma, *lead, step, secondaries);
                maxRaylErr = std::max(maxRaylErr, std::abs(gamma.getKineticEnergy() - e) / e + step.totalEnergyDeposit / e);
                raylSecondaries += static_cast<long>(secondaries.size());
                const Vector3D d = gamma.getDirection();
                maxDirErr = std::max(maxDirErr, std::abs(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z) - 1.0));
            }
            // Compton postStepDoIt
            {
                Track gamma(&gammaDef, Vector3D(), Vector3D(0.0, 0.0, 1.0), e);
                Step step;
                secondaries.clear();
                compt.postStepDoIt(gamma, *lead, step, secondaries);
                double sumE = gamma.getKineticEnergy() + step.totalEnergyDeposit;
                for (const Track& t : secondaries) {
                    sumE += t.getKineticEnergy();
                }
                maxComptErr = std::max(maxComptErr, std::abs(sumE - e) / e);
                const Vector3D d = gamma.getDirection();
                maxDirErr = std::max(maxDirErr, std::abs(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z) - 1.0));
            }
        }
    }

    if (maxRaylErr >= 1.0e-12 || maxComptErr >= 1.0e-12 || maxDirErr >= 1.0e-12 || raylSecondaries != 0) {
        ok = false;
    }

    std::printf("Conservation (%d events/energy): Rayleigh max|dE|/E = %.2e (sec = %ld), "
                "Compton max|dE|/E = %.2e, max|dir-1| = %.2e\n\n",
                nEvents, maxRaylErr, raylSecondaries, maxComptErr, maxDirErr);

    struct SlabCase {
        double eKeV;
        double thicknessCm;
    };
    const SlabCase cases[] = {{100.0, 0.02}, {300.0, 0.20}, {1000.0, 1.00}};

    std::printf("Thin lead slab transmission (%d photons each; phot + bound compt + Rayl):\n", nPhotons);
    std::printf("%8s %7s | %9s %9s %7s | %9s %9s %7s\n",
                "E0[keV]", "L [cm]", "T_MC", "exp(-uL)", "pull", "f_Rayl", "ref_Rayl", "pull");

    for (const auto& c : cases) {
        const double e0 = c.eKeV * Units::keV;
        const double muPhot = phot.computeMassAttenuation(e0, *lead) * lead->getDensity();
        const double muCompt = lead->getAtomDensity() * compt.computeCrossSectionPerAtom(e0, *lead);
        const double muRayl = lead->getAtomDensity() * rayl.computeCrossSectionPerAtom(e0, *lead);
        const double muTot = muPhot + muCompt + muRayl;

        const double expectedT = std::exp(-muTot * c.thicknessCm);
        const double expectedFirstRayl = (1.0 - expectedT) * (muRayl / muTot);

        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        const double hz = 0.5 * c.thicknessCm;
        rm.setWorldVolume(std::make_unique<Box>("World", 31.0, 31.0, hz + 1.0, vacuum));
        rm.addVolume(std::make_unique<Box>("LeadSlab", 30.0, 30.0, hz, lead));
        rm.addPhysicsProcess(std::make_unique<Ionization>());
        rm.addPhysicsProcess(std::make_unique<PhotoElectricEffect>(&electronDef));
        rm.addPhysicsProcess(std::make_unique<ComptonScattering>(&electronDef));
        rm.addPhysicsProcess(std::make_unique<RayleighScattering>());
        rm.setEnergyCut(10.0 * Units::keV);
        rm.setVerboseLevel(0);

        bool firstInteracted = false;
        long unscattered = 0;
        long firstRayl = 0;

        rm.setUserSteppingAction([&](const Track&, const Step& step) {
            if (!firstInteracted &&
                (step.processDefinedStep == "phot" || step.processDefinedStep == "compt" ||
                 step.processDefinedStep == "Rayl")) {
                firstInteracted = true;
                if (step.processDefinedStep == "Rayl") {
                    ++firstRayl;
                }
            }
        });

        for (int i = 0; i < nPhotons; ++i) {
            firstInteracted = false;
            Track photon(&gammaDef, Vector3D(0.0, 0.0, -(hz + 0.5)), Vector3D(0.0, 0.0, 1.0), e0);
            rm.simulateEvent(photon, 1.0 * Units::m, 10000);
            if (!firstInteracted) {
                ++unscattered;
            }
        }

        const double mcT = static_cast<double>(unscattered) / nPhotons;
        const double pullT = (mcT - expectedT) / std::sqrt(expectedT * (1.0 - expectedT) / nPhotons);

        const double mcFirstRayl = static_cast<double>(firstRayl) / nPhotons;
        const double sigmaRayl = std::sqrt(std::max(1.0e-12, expectedFirstRayl * (1.0 - expectedFirstRayl)) / nPhotons);
        const double pullRayl = (mcFirstRayl - expectedFirstRayl) / sigmaRayl;

        if (std::abs(pullT) >= kMaxPull || std::abs(pullRayl) >= kMaxPull) {
            ok = false;
        }
        std::printf("%8.0f %7.2f | %9.5f %9.5f %+7.2f | %9.5f %9.5f %+7.2f\n",
                    c.eKeV, c.thicknessCm, mcT, expectedT, pullT, mcFirstRayl, expectedFirstRayl, pullRayl);
    }

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
    allPassed &= testScatteringFactorsAndIntegrals();
    allPassed &= testCrossSectionsVsXcom();
    allPassed &= testRayleighSampling(quick ? 20000 : 200000);
    allPassed &= testBoundComptonSampling(quick ? 10000 : 100000);
    allPassed &= testConservationAndTransmission(quick ? 2000 : 10000, quick ? 10000 : 100000);

    std::printf("\n==============================================================================\n");
    std::printf("  OVERALL: %s\n", allPassed ? "PASS" : "FAIL");
    std::printf("==============================================================================\n");
    return allPassed ? 0 : 1;
}