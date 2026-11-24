// ==============================================================================
// [STAGE 8] Validation of stopped negative muon capture and bound decay-in-orbit.
//
// References:
//   - Capture rates, Z_eff, and Huff factors: Suzuki, Measday & Roalsvig,
//     Phys. Rev. C 35 (1987) 2212.
//   - Lead neutron multiplicity (<n> ~ 1.64-1.71): Measday, Phys. Rep. 354 (2001) 243.
//
//   TEST 1  Suzuki et al. capture rates, Huff factors, effective lifetimes & X-ray lines
//   TEST 2  Stopped mu- effective lifetime & capture fraction (Pb, C) and mu+ non-capture
//   TEST 3  Prompt K_alpha and L_alpha muonic X-ray yields on every stop (capture & DIO)
//   TEST 4  Capture final states: H two-body kinematics, Pb <n> & residual gamma bound
//   TEST 5  Event-by-event energy & DIO momentum conservation (with/without neutronDef)
//
// Run:  ./ValidationMuonCaptureApp
//       ./ValidationMuonCaptureApp quick
//       ./ValidationMuonCaptureApp 42
// ==============================================================================

#include "Box.h"
#include "Decay.h"
#include "Ionization.h"
#include "Material.h"
#include "MuonCapture.h"
#include "Particle.h"
#include "ParticleTable.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using PhysicalConstants::m_mu;

namespace {

constexpr double kMaxPull = 3.0;
constexpr double kFreeMuonLifetimeNs = 2196.9811;
constexpr double kProtonMassGeV = 0.93827208816;
constexpr double kNeutronMassGeV = 0.93956542052;

const Particle neutronDef(2112, "neutron", kNeutronMassGeV * Units::GeV, 0.0, 0.5, -1.0);

std::shared_ptr<Material> makeLead() {
    const DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

std::shared_ptr<Material> makeCarbon() {
    const DensityEffectData d{3.1550, 0.0480, 2.5387, 0.20762, 2.9532, 0.14};
    return std::make_shared<Material>("Carbon (C)", 2.00, 6.0, 12.011, 78.0 * Units::eV, 21.35 * Units::cm, d);
}

std::shared_ptr<Material> makeIron() {
    const DensityEffectData d{4.2911, -0.0012, 3.1531, 0.14680, 2.9632, 0.12};
    return std::make_shared<Material>("Iron (Fe)", 7.874, 26.0, 55.845, 286.0 * Units::eV, 1.757 * Units::cm, d);
}

std::shared_ptr<Material> makeHydrogen() {
    const DensityEffectData d{3.2632, 0.4759, 1.9215, 0.13483, 5.6249, 0.0};
    return std::make_shared<Material>("Hydrogen (H)", 0.0708, 1.0, 1.008, 21.8 * Units::eV, 890.4 * Units::cm, d);
}

void printHeader(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

// ==============================================================================
// TEST 1  Suzuki et al. (1987) capture rates, Huff factors & effective lifetimes
// ==============================================================================
bool testRatesAndConstants() {
    printHeader("TEST 1  Suzuki et al. (1987) capture rates, Huff factors & lifetimes");

    struct Ref {
        const char* name;
        int z;
        double lambda1e6;   // [10^6 s^-1]
        double refTauNs;    // [ns]
        double refCapProb;
    };

    const Ref refs[] = {
        {"C",  6,  0.0388, 2028.0, 0.0789},
        {"Al", 13, 0.7054,  864.0, 0.6095},
        {"Fe", 26, 4.411,   206.0, 0.9086},
        {"Pb", 82, 13.45,    72.3, 0.9722},
    };

    bool ok = true;
    double prevHuff = 2.0;

    std::printf("%-4s %4s %10s %10s %8s %9s %9s %8s %8s\n",
                "Elem", "Z", "Lam_c[us]", "Ref[us]", "Huff Q", "tau[ns]", "P_cap", "Ka[MeV]", "La[MeV]");

    for (const auto& r : refs) {
        const double lamUs = MuonCapture::getCaptureRatePerNs(r.z) * 1.0e3;
        const double q = MuonCapture::getHuffFactor(r.z);
        const double tauNs = MuonCapture::getEffectiveLifetimeNs(r.z);
        const double pCap = MuonCapture::getCaptureProbability(r.z);
        const double kAlphaMeV = MuonCapture::getKAlphaEnergy(r.z) / Units::MeV;
        const double lAlphaMeV = MuonCapture::getLAlphaEnergy(r.z) / Units::MeV;

        const double expectedTau = 1.0 / (MuonCapture::getCaptureRatePerNs(r.z) + q / kFreeMuonLifetimeNs);
        const double expectedProb = MuonCapture::getCaptureRatePerNs(r.z) * expectedTau;

        if (std::abs(lamUs / r.lambda1e6 - 1.0) > 1.0e-3 ||
            std::abs(tauNs - expectedTau) > 1.0e-9 ||
            std::abs(pCap - expectedProb) > 1.0e-12 ||
            std::abs(tauNs / r.refTauNs - 1.0) > 0.005 ||
            std::abs(pCap - r.refCapProb) > 0.002 ||
            q <= 0.0 || q > 1.0 || q >= prevHuff ||
            kAlphaMeV <= lAlphaMeV || lAlphaMeV <= 0.0) {
            ok = false;
        }
        prevHuff = q;

        std::printf("%-4s %4d %10.4f %10.4f %8.4f %9.2f %9.4f %8.3f %8.3f\n",
                    r.name, r.z, lamUs, r.lambda1e6, q, tauNs, pCap, kAlphaMeV, lAlphaMeV);
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 2  Stopped mu- lifetime & capture fraction (Pb, C) and mu+ non-capture
// ==============================================================================
bool testLifetimeAndBranching(const ParticleTable& table, int nEvents) {
    printHeader("TEST 2  Stopped mu- lifetime & capture fraction (Pb, C) and mu+ non-capture");

    const Particle& gammaDef = *table.findParticle("gamma");
    const Particle& electronDef = *table.findParticle("e-");
    const Particle& nuMuDef = *table.findParticle("nu_mu");
    const Particle& antiNuEDef = *table.findParticle("anti_nu_e");
    const Particle& muMinusDef = *table.findParticle("mu-");
    const Particle& muPlusDef = *table.findParticle("mu+");

    struct Case {
        const char* label;
        const Particle* particle;
        std::shared_ptr<Material> material;
        double expectedTauNs;
        double expectedCapProb;
    };

    const Case cases[] = {
        {"mu- in Pb", &muMinusDef, makeLead(),
         MuonCapture::getEffectiveLifetimeNs(82), MuonCapture::getCaptureProbability(82)},
        {"mu- in C",  &muMinusDef, makeCarbon(),
         MuonCapture::getEffectiveLifetimeNs(6),  MuonCapture::getCaptureProbability(6)},
        {"mu+ in Pb", &muPlusDef,  makeLead(),
         muPlusDef.getLifetime(), 0.0},
    };

    bool ok = true;
    std::printf("(%d stopped muons per case via full RunManager AtRest loop)\n", nEvents);
    std::printf("%-10s | %9s %9s %7s | %9s %9s %7s\n",
                "Case", "tau MC", "tau ref", "pull", "P_c MC", "P_c ref", "pull");

    for (const auto& c : cases) {
        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        rm.setWorldVolume(std::make_unique<Box>("World", 10.0 * Units::cm, 10.0 * Units::cm, 10.0 * Units::cm, vacuum));
        rm.addVolume(std::make_unique<Box>("Target", 5.0 * Units::cm, 5.0 * Units::cm, 5.0 * Units::cm, c.material));
        rm.addPhysicsProcess(std::make_unique<Ionization>());
        rm.addPhysicsProcess(std::make_unique<Decay>());
        rm.addPhysicsProcess(std::make_unique<MuonCapture>(&gammaDef, &electronDef, &nuMuDef, &antiNuEDef, &neutronDef));
        // Absorbs secondary e-/gamma immediately; never cuts mu- or mu+
        rm.setEnergyCut(1.0 * Units::GeV);
        rm.setVerboseLevel(0);

        bool stopped = false;
        double stopTime = 0.0;
        long nStopped = 0;
        long nCapture = 0;
        long nDecay = 0;
        double sumRestTime = 0.0;

        rm.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (std::abs(track.getParticleDefinition()->getPdgId()) != 13) {
                return;
            }
            if (!stopped && track.getTrackStatus() == TrackStatus::trackStopButAlive) {
                stopped = true;
                stopTime = track.getGlobalTime();
                return;
            }
            if (stopped && !track.isAlive()) {
                ++nStopped;
                sumRestTime += (track.getGlobalTime() - stopTime);
                if (step.processDefinedStep == "muMinusCapture") {
                    ++nCapture;
                } else if (step.processDefinedStep == "Decay") {
                    ++nDecay;
                }
            }
        });

        for (int i = 0; i < nEvents; ++i) {
            stopped = false;
            stopTime = 0.0;
            Track muon(c.particle, Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0), 1.0 * Units::MeV);
            rm.simulateEvent(muon, 1.0 * Units::cm, 1000);
        }

        if (nStopped == 0 || nCapture + nDecay != nStopped) {
            ok = false;
        }

        const double mcTau = sumRestTime / static_cast<double>(nStopped);
        const double tauPull = (mcTau - c.expectedTauNs) / (c.expectedTauNs / std::sqrt(static_cast<double>(nStopped)));

        const double mcProb = static_cast<double>(nCapture) / static_cast<double>(nStopped);
        double probPull = 0.0;
        if (c.expectedCapProb > 0.0 && c.expectedCapProb < 1.0) {
            const double sigmaP = std::sqrt(c.expectedCapProb * (1.0 - c.expectedCapProb) / static_cast<double>(nStopped));
            probPull = (mcProb - c.expectedCapProb) / sigmaP;
            if (std::abs(probPull) >= kMaxPull) {
                ok = false;
            }
        } else if (nCapture != 0) {
            ok = false;
        }

        if (std::abs(tauPull) >= kMaxPull) {
            ok = false;
        }

        std::printf("%-10s | %9.2f %9.2f %+7.2f | %9.5f %9.5f %+7.2f\n",
                    c.label, mcTau, c.expectedTauNs, tauPull, mcProb, c.expectedCapProb, probPull);
    }

    const double noHuffPb = MuonCapture::getCaptureRatePerNs(82) /
                            (MuonCapture::getCaptureRatePerNs(82) + 1.0 / kFreeMuonLifetimeNs);
    std::printf("Note: Without Huff factor Q(Pb)=%.4f, Pb P_c would be %.5f (vs %.5f with Huff).\n",
                MuonCapture::getHuffFactor(82), noHuffPb, MuonCapture::getCaptureProbability(82));
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 3  Prompt K_alpha and L_alpha muonic X-ray yields on every stop
// ==============================================================================
bool testMuonicXrays(const ParticleTable& table, int nEvents) {
    printHeader("TEST 3  Prompt K_alpha & L_alpha X-ray yields on EVERY stop (capture & DIO)");

    const Particle& gammaDef = *table.findParticle("gamma");
    const Particle& electronDef = *table.findParticle("e-");
    const Particle& nuMuDef = *table.findParticle("nu_mu");
    const Particle& antiNuEDef = *table.findParticle("anti_nu_e");
    const Particle& muMinusDef = *table.findParticle("mu-");
    const Particle& muPlusDef = *table.findParticle("mu+");

    struct Case {
        const char* label;
        const Particle* particle;
        std::shared_ptr<Material> material;
        int z;
        bool expectXrays;
    };

    const Case cases[] = {
        {"mu- in C",  &muMinusDef, makeCarbon(), 6,  true},
        {"mu- in Pb", &muMinusDef, makeLead(),   82, true},
        {"mu+ in Pb", &muPlusDef,  makeLead(),   82, false},
    };

    bool ok = true;
    std::printf("%-10s | %7s %8s %8s | %7s %8s %8s\n",
                "Case", "N_cap", "Ka/cap", "La/cap", "N_dio", "Ka/dio", "La/dio");

    for (const auto& c : cases) {
        const double eKa = MuonCapture::getKAlphaEnergy(c.z);
        const double eLa = MuonCapture::getLAlphaEnergy(c.z);
        const double tol = 1.0 * Units::eV;

        RunManager rm;
        auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
        rm.setWorldVolume(std::make_unique<Box>("World", 10.0 * Units::cm, 10.0 * Units::cm, 10.0 * Units::cm, vacuum));
        rm.addVolume(std::make_unique<Box>("Target", 5.0 * Units::cm, 5.0 * Units::cm, 5.0 * Units::cm, c.material));
        rm.addPhysicsProcess(std::make_unique<Ionization>());
        rm.addPhysicsProcess(std::make_unique<Decay>());
        rm.addPhysicsProcess(std::make_unique<MuonCapture>(&gammaDef, &electronDef, &nuMuDef, &antiNuEDef, &neutronDef));
        rm.setEnergyCut(1.0 * Units::GeV);
        rm.setVerboseLevel(0);

        bool stopped = false;
        long nCap = 0, kaCap = 0, laCap = 0;
        long nDio = 0, kaDio = 0, laDio = 0;

        rm.setUserSteppingAction([&](const Track& track, const Step& step) {
            if (std::abs(track.getParticleDefinition()->getPdgId()) != 13) {
                return;
            }
            if (!stopped && track.getTrackStatus() == TrackStatus::trackStopButAlive) {
                stopped = true;
                return;
            }
            if (!stopped || track.isAlive() || step.secondaryInCurrentStep == nullptr) {
                return;
            }

            long ka = 0, la = 0;
            for (const Track& sec : *step.secondaryInCurrentStep) {
                if (sec.getParticleDefinition()->getPdgId() != 22) {
                    continue;
                }
                const double e = sec.getKineticEnergy();
                if (std::abs(e - eKa) < tol) ++ka;
                if (std::abs(e - eLa) < tol) ++la;
            }

            if (step.processDefinedStep == "muMinusCapture") {
                ++nCap;
                kaCap += ka;
                laCap += la;
            } else if (step.processDefinedStep == "Decay") {
                ++nDio;
                kaDio += ka;
                laDio += la;
            }
        });

        for (int i = 0; i < nEvents; ++i) {
            stopped = false;
            Track muon(c.particle, Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0), 1.0 * Units::MeV);
            rm.simulateEvent(muon, 1.0 * Units::cm, 1000);
        }

        const double kaCapRate = (nCap > 0) ? static_cast<double>(kaCap) / static_cast<double>(nCap) : 0.0;
        const double laCapRate = (nCap > 0) ? static_cast<double>(laCap) / static_cast<double>(nCap) : 0.0;
        const double kaDioRate = (nDio > 0) ? static_cast<double>(kaDio) / static_cast<double>(nDio) : 0.0;
        const double laDioRate = (nDio > 0) ? static_cast<double>(laDio) / static_cast<double>(nDio) : 0.0;

        if (c.expectXrays) {
            if (nCap == 0 || nDio == 0 || kaCap != nCap || laCap != nCap || kaDio != nDio || laDio != nDio) {
                ok = false;
            }
        } else {
            if (nCap != 0 || kaDio != 0 || laDio != 0) {
                ok = false;
            }
        }

        std::printf("%-10s | %7ld %8.4f %8.4f | %7ld %8.4f %8.4f\n",
                    c.label, nCap, kaCapRate, laCapRate, nDio, kaDioRate, laDioRate);
    }

    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 4  Capture final states: H 2-body kinematics, Pb <n> & residual gamma bound
// ==============================================================================
bool testCaptureSpectraAndMultiplicity(const ParticleTable& table, int nEvents) {
    printHeader("TEST 4  Capture final states: H 2-body, Pb neutron multiplicity & gamma bound");

    const Particle& gammaDef = *table.findParticle("gamma");
    const Particle& electronDef = *table.findParticle("e-");
    const Particle& nuMuDef = *table.findParticle("nu_mu");
    const Particle& antiNuEDef = *table.findParticle("anti_nu_e");

    // Long-lived mu- probe ensures nuclear capture dominates even in Z = 1 (kinematics use m_mu)
    const Particle captureProbe(13, "mu-_probe", m_mu, -1.0, 0.5, 1.0e12);

    bool ok = true;

    // 4a. Hydrogen exact 2-body: mu- + p -> n + nu_mu
    {
        auto hydrogen = makeHydrogen();
        MuonCapture cap(&gammaDef, &electronDef, &nuMuDef, &antiNuEDef, &neutronDef);
        cap.setVerboseLevel(0);

        const double w = kProtonMassGeV * Units::GeV + m_mu - MuonCapture::get1sBindingEnergy(1);
        const double mn = kNeutronMassGeV * Units::GeV;
        const double refEnu = 0.5 * (w * w - mn * mn) / w;
        const double refTn = w - mn - refEnu;

        double maxNuErr = 0.0;
        double maxNErr = 0.0;
        double maxBackToBackErr = 0.0;
        long nCap = 0;

        std::vector<Track> sec;
        sec.reserve(8);
        for (int i = 0; i < 2000; ++i) {
            Track mu(&captureProbe, Vector3D(), Vector3D(0.0, 0.0, 1.0), 0.0);
            Step step;
            sec.clear();
            cap.atRestDoIt(mu, *hydrogen, step, sec);
            if (step.processDefinedStep != "muMinusCapture") {
                continue;
            }
            ++nCap;
            double eNu = 0.0, tN = 0.0;
            Vector3D dNu, dN;
            for (const Track& t : sec) {
                const int pdg = t.getParticleDefinition()->getPdgId();
                if (pdg == 14) {
                    eNu = t.getKineticEnergy();
                    dNu = t.getDirection();
                } else if (pdg == 2112) {
                    tN = t.getKineticEnergy();
                    dN = t.getDirection();
                }
            }
            maxNuErr = std::max(maxNuErr, std::abs(eNu - refEnu) / refEnu);
            maxNErr = std::max(maxNErr, std::abs(tN - refTn) / refTn);
            const double dot = dNu.x * dN.x + dNu.y * dN.y + dNu.z * dN.z;
            maxBackToBackErr = std::max(maxBackToBackErr, std::abs(dot + 1.0));
        }

        if (nCap == 0 || maxNuErr > 1.0e-6 || maxNErr > 1.0e-6 || maxBackToBackErr > 1.0e-12) {
            ok = false;
        }
        std::printf("H  2-body (mu- p -> n nu_mu): E_nu = %.3f MeV (err %.1e), T_n = %.3f MeV (err %.1e), |u_nu.u_n + 1| = %.1e\n",
                    refEnu / Units::MeV, maxNuErr, refTn / Units::MeV, maxNErr, maxBackToBackErr);
    }

    // 4b. Lead neutron multiplicity and de-excitation gamma cap (with & without neutronDef)
    {
        auto lead = makeLead();
        const double eKa = MuonCapture::getKAlphaEnergy(82);
        const double eLa = MuonCapture::getLAlphaEnergy(82);
        const double tol = 1.0 * Units::eV;
        constexpr double kMaxAllowedPbResidualGammaMeV = 6.846 + 1.0e-3; // max(S_n1, S_n2) for 208Tl -> 207Tl

        double meanNWithDef = 0.0;
        double varNWithDef = 0.0;

        std::printf("\n%-16s %8s %9s %10s %11s %12s\n",
                    "Pb Mode", "<n>", "<T_n>MeV", "<E_nu>MeV", "nTracks/cap", "maxE_gam MeV");

        for (bool hasNeutronDef : {true, false}) {
            MuonCapture cap(&gammaDef, &electronDef, &nuMuDef, &antiNuEDef, hasNeutronDef ? &neutronDef : nullptr);
            cap.setVerboseLevel(0);

            long nCap = 0;
            long sumN = 0;
            long sumN2 = 0;
            long sumNeutronTracks = 0;
            double sumTn = 0.0;
            double sumEnu = 0.0;
            double maxResidualGamma = 0.0;

            std::vector<Track> sec;
            sec.reserve(12);

            for (int i = 0; i < nEvents; ++i) {
                Track mu(&captureProbe, Vector3D(), Vector3D(0.0, 0.0, 1.0), 0.0);
                Step step;
                sec.clear();
                cap.atRestDoIt(mu, *lead, step, sec);
                if (step.processDefinedStep != "muMinusCapture") {
                    continue;
                }
                ++nCap;
                const int mult = cap.getLastNeutronCount();
                sumN += mult;
                sumN2 += static_cast<long>(mult) * mult;

                for (const Track& t : sec) {
                    const int pdg = t.getParticleDefinition()->getPdgId();
                    const double e = t.getKineticEnergy();
                    if (pdg == 2112) {
                        ++sumNeutronTracks;
                        sumTn += e;
                    } else if (pdg == 14) {
                        sumEnu += e;
                    } else if (pdg == 22) {
                        if (std::abs(e - eKa) > tol && std::abs(e - eLa) > tol) {
                            maxResidualGamma = std::max(maxResidualGamma, e);
                        }
                    }
                }
            }

            const double meanN = static_cast<double>(sumN) / static_cast<double>(nCap);
            const double varN = std::max(0.0, static_cast<double>(sumN2) / static_cast<double>(nCap) - meanN * meanN);
            const double meanTnMeV = (sumNeutronTracks > 0) ? (sumTn / static_cast<double>(sumNeutronTracks)) / Units::MeV : 0.0;
            const double meanEnuMeV = (sumEnu / static_cast<double>(nCap)) / Units::MeV;
            const double nTracksPerCap = static_cast<double>(sumNeutronTracks) / static_cast<double>(nCap);
            const double maxGamMeV = maxResidualGamma / Units::MeV;

            if (hasNeutronDef) {
                meanNWithDef = meanN;
                varNWithDef = varN;
                // Measured natural Pb multiplicity is ~1.64-1.71 neutrons/capture
                if (meanN < 1.45 || meanN > 1.95 || sumNeutronTracks != sumN) {
                    ok = false;
                }
            } else {
                const double sigmaDiff = std::sqrt((varNWithDef + varN) / static_cast<double>(nCap));
                const double pull = (sigmaDiff > 0.0) ? (meanN - meanNWithDef) / sigmaDiff : 0.0;
                if (sumNeutronTracks != 0 || std::abs(pull) >= kMaxPull) {
                    ok = false;
                }
            }

            if (maxGamMeV > kMaxAllowedPbResidualGammaMeV) {
                ok = false;
            }

            std::printf("%-16s %8.3f %9.3f %10.2f %11.3f %12.3f\n",
                        hasNeutronDef ? "neutronDef ON" : "neutronDef NULL",
                        meanN, meanTnMeV, meanEnuMeV, nTracksPerCap, maxGamMeV);
        }
    }

    std::printf("Pb <n> reference: ~1.64-1.71 n/capture; residual nuclear gamma must stay <= S_n2 (6.846 MeV).\n");
    std::printf("Result: %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// ==============================================================================
// TEST 5  Event-by-event energy & DIO momentum conservation
// ==============================================================================
bool testEnergyAndMomentumConservation(const ParticleTable& table, int nEvents) {
    printHeader("TEST 5  Event-by-event energy & DIO momentum conservation");

    const Particle& gammaDef = *table.findParticle("gamma");
    const Particle& electronDef = *table.findParticle("e-");
    const Particle& nuMuDef = *table.findParticle("nu_mu");
    const Particle& antiNuEDef = *table.findParticle("anti_nu_e");
    const Particle& muMinusDef = *table.findParticle("mu-");
    const Particle captureProbe(13, "mu-_probe", m_mu, -1.0, 0.5, 1.0e12);

    const std::vector<std::shared_ptr<Material>> materials = {
        makeHydrogen(), makeCarbon(), makeIron(), makeLead()
    };

    struct Mode {
        const char* label;
        bool hasNeutron;
        bool emitXrays;
    };
    const Mode modes[] = {
        {"n=ON,  Xray=ON",  true,  true},
        {"n=OFF, Xray=ON",  false, true},
        {"n=ON,  Xray=OFF", true,  false},
    };

    double overallMaxCapErr = 0.0;
    double overallMaxDioErr = 0.0;
    double overallMaxMomErr = 0.0;
    double overallMinLocal = 1.0e9;

    std::printf("%-14s %-16s | %11s %11s %11s %12s\n",
                "Material", "Mode", "cap |dE|/m", "dio |dE|/m", "dio |dp|/m", "min local");

    std::vector<Track> sec;
    sec.reserve(12);

    for (const auto& mat : materials) {
        for (const auto& mode : modes) {
            MuonCapture cap(&gammaDef, &electronDef, &nuMuDef, &antiNuEDef,
                            mode.hasNeutron ? &neutronDef : nullptr);
            cap.setMuonicXrays(mode.emitXrays);
            cap.setVerboseLevel(0);

            double maxCapErr = 0.0;
            double maxDioErr = 0.0;
            double maxMomErr = 0.0;
            double minLocal = 1.0e9;

            // Half real mu- (exercises DIO + capture), half long-lived probe (ensures capture even in H/C)
            for (int pass = 0; pass < 2; ++pass) {
                const Particle* pDef = (pass == 0) ? &muMinusDef : &captureProbe;
                const int nPass = nEvents / 2;

                for (int i = 0; i < nPass; ++i) {
                    Track mu(pDef, Vector3D(), Vector3D(0.0, 0.0, 1.0), 0.0);
                    Step step;
                    sec.clear();
                    cap.atRestDoIt(mu, *mat, step, sec);

                    minLocal = std::min(minLocal, step.totalEnergyDeposit);

                    if (step.processDefinedStep == "muMinusCapture") {
                        double sumE = step.totalEnergyDeposit + cap.getLastNuclearMassChange();
                        for (const Track& t : sec) {
                            sumE += t.getKineticEnergy();
                        }
                        maxCapErr = std::max(maxCapErr, std::abs(sumE - m_mu) / m_mu);
                    } else if (step.processDefinedStep == "Decay") {
                        double sumE = step.totalEnergyDeposit;
                        double px = 0.0, py = 0.0, pz = 0.0;
                        for (const Track& t : sec) {
                            const double eKin = t.getKineticEnergy();
                            const double mass = t.getMass();
                            sumE += eKin + mass;
                            if (t.getParticleDefinition()->getPdgId() != 22) {
                                const double pMag = std::sqrt(std::max(0.0, eKin * (eKin + 2.0 * mass)));
                                const Vector3D d = t.getDirection();
                                px += pMag * d.x;
                                py += pMag * d.y;
                                pz += pMag * d.z;
                            }
                        }
                        maxDioErr = std::max(maxDioErr, std::abs(sumE - m_mu) / m_mu);
                        const double pRes = std::sqrt(px * px + py * py + pz * pz);
                        maxMomErr = std::max(maxMomErr, pRes / m_mu);
                    }
                }
            }

            overallMaxCapErr = std::max(overallMaxCapErr, maxCapErr);
            overallMaxDioErr = std::max(overallMaxDioErr, maxDioErr);
            overallMaxMomErr = std::max(overallMaxMomErr, maxMomErr);
            overallMinLocal = std::min(overallMinLocal, minLocal);

            std::printf("%-14s %-16s | %11.2e %11.2e %11.2e %9.3f MeV\n",
                        mat->getName().c_str(), mode.label, maxCapErr, maxDioErr, maxMomErr, minLocal / Units::MeV);
        }
    }

    const bool ok = (overallMaxCapErr < 1.0e-12) &&
                    (overallMaxDioErr < 1.0e-12) &&
                    (overallMaxMomErr < 1.0e-12) &&
                    (overallMinLocal >= 0.0);

    if (overallMinLocal >= 0.0) {
        std::printf("min local >= 0: no event produces a negative local energy deposit.\n");
    } else {
        std::printf("FAIL: negative local energy deposit encountered (%.6f MeV).\n", overallMinLocal / Units::MeV);
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

    ParticleTable particleTable;

    bool allPassed = true;
    allPassed &= testRatesAndConstants();
    allPassed &= testLifetimeAndBranching(particleTable, quick ? 20000 : 100000);
    allPassed &= testMuonicXrays(particleTable, quick ? 5000 : 20000);
    allPassed &= testCaptureSpectraAndMultiplicity(particleTable, quick ? 5000 : 25000);
    allPassed &= testEnergyAndMomentumConservation(particleTable, quick ? 2000 : 10000);

    std::printf("\n==============================================================================\n");
    std::printf("  OVERALL: %s\n", allPassed ? "PASS" : "FAIL");
    std::printf("==============================================================================\n");
    return allPassed ? 0 : 1;
}