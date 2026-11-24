#include "MuonCapture.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using PhysicalConstants::m_e;
using PhysicalConstants::m_mu;

namespace {

constexpr double kFreeMuonLifetimeNs = 2196.9811; // PDG [ns]
constexpr double kFreeMuonDecayRatePerNs = 1.0 / kFreeMuonLifetimeNs;
constexpr double kAtomicMassUnitGeV = 0.93149410242;
constexpr double kProtonMassGeV = 0.93827208816;
constexpr double kNeutronMassGeV = 0.93956542052;

// Dominant stable isotope (Z, A) -> capture daughter (Z-1, A)
struct ElementCaptureEntry {
    int z;
    double zEff;             // Ford-Wills / Suzuki effective Z
    double huffQ;            // Tabulated Huff factor Q (Suzuki et al. 1987)
    double lambda1e6PerSec;  // Total capture rate [10^6 s^-1] (Suzuki et al. 1987)
    double binding1sMeV;     // Muonic 1s binding energy B_1s [MeV]
    double kAlphaMeV;        // 2p -> 1s muonic X-ray energy [MeV]
    double lAlphaMeV;        // 3d -> 2p muonic X-ray energy [MeV]
    double deltaMassNucMeV;  // M_nuc(Z-1, A) - M_nuc(Z, A) = Q_beta-(Z-1, A) + m_e [MeV]
    double sepN1MeV;         // 1st neutron separation energy S_n(Z-1, A) [MeV]
    double sepN2MeV;         // 2nd+ neutron separation energy S_n(Z-1, A-1) [MeV]
};

constexpr std::array<ElementCaptureEntry, 20> kCaptureTable = {{
    // Z,  zEff,   huffQ, Lambda_c,  B_1s,    K_a,     L_a,     dM_nuc,  S_n1,   S_n2    (Parent -> Daughter)
    { 1,  1.00,  1.000,   0.00046,  0.00253, 0.00190, 0.00035,  1.293,  0.000,  0.000}, // 1H   -> n
    { 6,  5.72,  0.998,   0.0388,   0.100,   0.075,   0.014,   13.880,  3.370, 11.454}, // 12C  -> 12B
    { 7,  6.61,  0.997,   0.0693,   0.137,   0.102,   0.019,    0.667,  8.176,  4.946}, // 14N  -> 14C
    { 8,  7.49,  0.997,   0.1026,   0.178,   0.133,   0.025,   10.932,  2.489, 10.833}, // 16O  -> 16N
    {13, 11.48,  0.993,   0.7054,   0.463,   0.347,   0.063,    3.121,  6.443, 11.093}, // 27Al -> 27Mg
    {14, 12.22,  0.992,   0.8712,   0.535,   0.400,   0.073,    5.153,  7.725, 13.058}, // 28Si -> 28Al
    {18, 14.89,  0.988,   1.41,     0.862,   0.644,   0.119,    7.990,  5.820,  8.070}, // 40Ar -> 40Cl
    {20, 16.15,  0.985,   2.557,    1.055,   0.784,   0.146,    1.822,  7.800, 13.077}, // 40Ca -> 40K
    {26, 19.59,  0.975,   4.411,    1.720,   1.255,   0.254,    4.206,  7.271, 10.226}, // 56Fe -> 56Mn
    {29, 21.12,  0.968,   5.676,    2.090,   1.512,   0.315,    0.578,  6.838, 10.596}, // 63Cu -> 63Ni
    {32, 22.43,  0.960,   5.57,     2.480,   1.775,   0.380,    5.880,  6.420,  9.190}, // 74Ge -> 74Ga
    {47, 27.95,  0.925,  10.86,     4.730,   3.170,   0.825,    0.545,  6.539,  9.561}, // 107Ag-> 107Pd
    {50, 28.64,  0.916,  10.44,     5.210,   3.440,   0.935,    5.880,  5.300,  8.150}, // 120Sn-> 120In
    {53, 29.27,  0.910,  11.20,     5.690,   3.710,   1.050,    1.209,  6.288,  9.114}, // 127I -> 127Te
    {54, 29.46,  0.908,  11.10,     5.860,   3.800,   1.090,    4.090,  6.450,  8.580}, // 132Xe-> 132I
    {74, 33.03,  0.861,  12.36,     9.220,   5.330,   2.100,    3.380,  5.650,  7.230}, // 184W -> 184Ta
    {79, 33.64,  0.850,  13.07,    10.080,   5.710,   2.380,    1.230,  5.846,  7.922}, // 197Au-> 197Pt
    {82, 34.18,  0.844,  13.45,    10.590,   5.900,   2.550,    5.510,  3.787,  6.846}, // 208Pb-> 208Tl
    {83, 34.00,  0.842,  13.26,    10.760,   6.020,   2.590,    1.155,  3.937,  7.368}, // 209Bi-> 209Pb
    {92, 34.94,  0.818,  12.59,    12.200,   6.450,   3.170,    4.470,  4.500,  6.200}, // 238U -> 238Pa
}};

// Exact Z lookup: nuclear-structure values (DeltaM, S_n) cannot be interpolated across nuclei
const ElementCaptureEntry& lookupElement(int z) {
    for (const auto& entry : kCaptureTable) {
        if (entry.z == z) {
            return entry;
        }
    }
    throw std::runtime_error("MuonCapture: element Z = " + std::to_string(z) + " is not in kCaptureTable");
}

// Exact, unbiased sampler for P(T) ~ T * exp(-T / theta) on [0, tMax].
// Acceptance >= 39% for all tMax > 0 (approaches 100% at both extremes).
double sampleMaxwellian(double theta, double tMax) {
    if (theta <= 0.0 || tMax <= 0.0) {
        return 0.0;
    }
    if (tMax <= 1.5 * theta) {
        for (;;) {
            const double t = tMax * std::sqrt(RandomEngine::uniNumGen());
            if (RandomEngine::uniNumGen() <= std::exp(-t / theta)) {
                return t;
            }
        }
    }
    for (;;) {
        const double u1 = 1.0 - RandomEngine::uniNumGen();
        const double u2 = 1.0 - RandomEngine::uniNumGen();
        const double t = -theta * std::log(u1 * u2);
        if (t <= tMax) {
            return t;
        }
    }
}

} // namespace

// ==============================================================================
// Static Element Physics Lookups
// ==============================================================================
Vector3D MuonCapture::sampleIsotropicDirection() {
    const double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    return Vector3D(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

double MuonCapture::getFreeDecayRatePerNs(const Track& track) {
    const double tau = track.getParticleDefinition()->getLifetime();
    return (tau > 0.0) ? (1.0 / tau) : kFreeMuonDecayRatePerNs;
}

double MuonCapture::getEffectiveZ(int z) {
    return (z >= 1) ? lookupElement(z).zEff : 0.0;
}

double MuonCapture::getHuffFactor(int z) {
    return (z >= 1) ? lookupElement(z).huffQ : 1.0;
}

double MuonCapture::getCaptureRatePerNs(int z) {
    return (z >= 1) ? lookupElement(z).lambda1e6PerSec * 1.0e-3 : 0.0;
}

double MuonCapture::getCaptureProbability(int z) {
    const double lambdaC = getCaptureRatePerNs(z);
    const double lambdaD = getHuffFactor(z) * kFreeMuonDecayRatePerNs;
    return (lambdaC > 0.0) ? lambdaC / (lambdaC + lambdaD) : 0.0;
}

double MuonCapture::getEffectiveLifetimeNs(int z) {
    const double totalRate = getCaptureRatePerNs(z) + getHuffFactor(z) * kFreeMuonDecayRatePerNs;
    return (totalRate > 0.0) ? 1.0 / totalRate : GeoDefs::kInfinity;
}

double MuonCapture::get1sBindingEnergy(int z) {
    return (z >= 1) ? lookupElement(z).binding1sMeV * Units::MeV : 0.0;
}

double MuonCapture::getKAlphaEnergy(int z) {
    return (z >= 1) ? lookupElement(z).kAlphaMeV * Units::MeV : 0.0;
}

double MuonCapture::getLAlphaEnergy(int z) {
    return (z >= 1) ? lookupElement(z).lAlphaMeV * Units::MeV : 0.0;
}

double MuonCapture::getNuclearMassChange(int z) {
    return (z >= 1) ? lookupElement(z).deltaMassNucMeV * Units::MeV : 0.0;
}

double MuonCapture::getFirstNeutronSeparationEnergy(int z) {
    return (z >= 1) ? lookupElement(z).sepN1MeV * Units::MeV : 0.0;
}

// ==============================================================================
// At-Rest GPIL: Total Stopped Lifetime (Capture + Bound DIO) [ns]
// ==============================================================================
double MuonCapture::atRestGetPhysicalInteractionLength(const Track& track, const Material& material) {
    const int z = static_cast<int>(std::lround(material.getZ()));
    if (z < 1) {
        return GeoDefs::kInfinity;
    }
    const double totalRate = getCaptureRatePerNs(z) + getHuffFactor(z) * getFreeDecayRatePerNs(track);
    return (totalRate > 0.0) ? -std::log(1.0 - RandomEngine::uniNumGen()) / totalRate : GeoDefs::kInfinity;
}

// ==============================================================================
// Muonic Cascade (emitted on EVERY stop, before capture or DIO)
// ==============================================================================
void MuonCapture::emitAtomicCascade(int z, const Vector3D& pos, Step& step,
                                    std::vector<Track>& secondaries) const {
    const ElementCaptureEntry& elem = lookupElement(z);
    const double binding1s = elem.binding1sMeV * Units::MeV;
    if (binding1s <= 0.0) {
        return;
    }

    if (!emitMuonicXrays || gammaDef == nullptr) {
        step.addTotalEnergyDeposit(binding1s);
        return;
    }

    const double kAlpha = std::min(binding1s, elem.kAlphaMeV * Units::MeV);
    const double lAlpha = std::min(binding1s - kAlpha, elem.lAlphaMeV * Units::MeV);
    const double outerCascade = std::max(0.0, binding1s - kAlpha - lAlpha);

    if (kAlpha > 0.0) {
        secondaries.emplace_back(gammaDef, pos, sampleIsotropicDirection(), kAlpha);
    }
    if (lAlpha > 0.0) {
        secondaries.emplace_back(gammaDef, pos, sampleIsotropicDirection(), lAlpha);
    }
    if (outerCascade > 0.0) {
        step.addTotalEnergyDeposit(outerCascade);
    }
}

// ==============================================================================
// Bound Decay-in-Orbit: mu-(1s) -> e- + anti_nu_e + nu_mu
// Simplification: samples the free Michel spectrum with M_eff = m_mu - B_1s,
// omitting nuclear-recoil smearing and the high-energy DIO endpoint tail.
// ==============================================================================
void MuonCapture::performBoundDecay(double binding1s, const Vector3D& pos, Step& step,
                                    std::vector<Track>& secondaries) {
    lastNuclearMassChange = 0.0;
    lastNeutronCount = 0;
    step.processDefinedStep = "Decay";

    const double mEff = m_mu - binding1s;
    if (mEff <= m_e || electronDef == nullptr) {
        step.addTotalEnergyDeposit(std::max(0.0, mEff));
        return;
    }

    const double eMax = 0.5 * (mEff * mEff + m_e * m_e) / mEff;
    const double xMin = m_e / eMax;
    double x = 1.0;
    do {
        x = std::cbrt(RandomEngine::uniNumGen());
    } while (x < xMin || RandomEngine::uniNumGen() * 3.0 > (3.0 - 2.0 * x));

    const double eElec = x * eMax;
    const double tElec = std::max(0.0, eElec - m_e);
    const double pElec = std::sqrt(std::max(0.0, eElec * eElec - m_e * m_e));
    const Vector3D eDir = sampleIsotropicDirection();
    secondaries.emplace_back(electronDef, pos, eDir, tElec);

    // Exact 2-neutrino decay of recoil system (E_vv = mEff - eElec, p_vv = -p_e)
    const double eNuNu = std::max(0.0, mEff - eElec);
    const double mNuNu = std::sqrt(std::max(0.0, eNuNu * eNuNu - pElec * pElec));
    const double cosStar = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinStar = std::sqrt(std::max(0.0, (1.0 - cosStar) * (1.0 + cosStar)));
    const double phiStar = 2.0 * M_PI * RandomEngine::uniNumGen();

    const double eNu1 = std::max(0.0, 0.5 * (eNuNu + pElec * cosStar));
    const double eNu2 = std::max(0.0, eNuNu - eNu1);
    const double pPar1 = 0.5 * (eNuNu * cosStar + pElec);
    const double pPerp = 0.5 * mNuNu * sinStar;

    const Vector3D oppDir = eDir * -1.0;
    Vector3D p1(pPerp * std::cos(phiStar), pPerp * std::sin(phiStar), pPar1);
    Vector3D p2(-pPerp * std::cos(phiStar), -pPerp * std::sin(phiStar), pElec - pPar1);
    p1.rotateUz(oppDir);
    p2.rotateUz(oppDir);

    if (nuMuDef != nullptr && eNu1 > 0.0) {
        secondaries.emplace_back(nuMuDef, pos, p1.unit(), eNu1);
    }
    if (antiNuEDef != nullptr && eNu2 > 0.0) {
        secondaries.emplace_back(antiNuEDef, pos, p2.unit(), eNu2);
    }

    if (verboseLevel > 0) {
        std::cout << "    [Decay/DIO] Z-bound | T_e=" << tElec / Units::MeV
                  << " MeV | E_vv=" << eNuNu / Units::MeV << " MeV\n";
    }
}

// ==============================================================================
// Nuclear Capture: mu-(1s) + (Z, A) -> (Z-1, A)* + nu_mu
// ==============================================================================
void MuonCapture::performNuclearCapture(int z, const Material& material, const Vector3D& pos,
                                        Step& step, std::vector<Track>& secondaries) {
    const ElementCaptureEntry& elem = lookupElement(z);
    const double binding1s = elem.binding1sMeV * Units::MeV;
    const double deltaM = elem.deltaMassNucMeV * Units::MeV;

    lastNuclearMassChange = deltaM;
    lastNeutronCount = 0;
    step.processDefinedStep = processName;

    // 1. Hydrogen: exact 2-body mu- + p -> n + nu_mu
    if (z == 1) {
        const double mp = kProtonMassGeV * Units::GeV;
        const double mn = kNeutronMassGeV * Units::GeV;
        lastNuclearMassChange = mn - mp;

        const double w = mp + m_mu - binding1s;
        const double eNu = std::max(0.0, 0.5 * (w * w - mn * mn) / w);
        const double tNeutron = std::max(0.0, w - mn - eNu);
        const Vector3D nuDir = sampleIsotropicDirection();

        if (nuMuDef != nullptr && eNu > 0.0) {
            secondaries.emplace_back(nuMuDef, pos, nuDir, eNu);
        }
        if (neutronDef != nullptr && tNeutron > 0.0) {
            secondaries.emplace_back(neutronDef, pos, nuDir * -1.0, tNeutron);
        } else if (tNeutron > 0.0) {
            step.addTotalEnergyDeposit(tNeutron);
        }
        lastNeutronCount = 1;
        return;
    }

    // 2. Complex nucleus: exact 2-body kinematics with M(Z-1, A) = M(Z, A) + deltaM
    const double eAvail = std::max(0.0, m_mu - binding1s - deltaM);
    const double nucMass = material.getA() * kAtomicMassUnitGeV * Units::GeV;
    constexpr double kExcitationScaleMeV = 7.5;
    const double eStar = sampleMaxwellian(kExcitationScaleMeV * Units::MeV, eAvail);

    const double w = nucMass + m_mu - binding1s;
    const double mStar = nucMass + deltaM + eStar;
    const double eNu = std::max(0.0, 0.5 * (w * w - mStar * mStar) / w);
    const double tRecoil = std::max(0.0, eAvail - eStar - eNu);

    if (nuMuDef != nullptr && eNu > 0.0) {
        secondaries.emplace_back(nuMuDef, pos, sampleIsotropicDirection(), eNu);
    }
    if (tRecoil > 0.0) {
        step.addTotalEnergyDeposit(tRecoil);
    }

    // 3. Evaporation neutrons (always runs, even when neutronDef == nullptr)
    double remainingExcitation = eStar;
    const double sn1 = elem.sepN1MeV * Units::MeV;
    const double sn2 = elem.sepN2MeV * Units::MeV;
    constexpr double kNuclearTempMeV = 1.4;

    while (true) {
        const double sn = (lastNeutronCount == 0) ? sn1 : sn2;
        if (sn <= 0.0 || remainingExcitation <= sn) {
            break;
        }
        const double tNeutron = sampleMaxwellian(kNuclearTempMeV * Units::MeV, remainingExcitation - sn);
        if (neutronDef != nullptr) {
            secondaries.emplace_back(neutronDef, pos, sampleIsotropicDirection(), tNeutron);
        } else {
            step.addTotalEnergyDeposit(tNeutron);
        }
        remainingExcitation -= (tNeutron + sn);
        lastNuclearMassChange += sn;
        ++lastNeutronCount;
    }

    // 4. Residual nuclear de-excitation gamma (guaranteed <= S_n)
    if (remainingExcitation > 0.0) {
        if (gammaDef != nullptr) {
            secondaries.emplace_back(gammaDef, pos, sampleIsotropicDirection(), remainingExcitation);
        } else {
            step.addTotalEnergyDeposit(remainingExcitation);
        }
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Z=" << z << " | E_nu=" << eNu / Units::MeV
                  << " MeV | E*=" << eStar / Units::MeV << " MeV | n=" << lastNeutronCount
                  << " | E_gamma=" << remainingExcitation / Units::MeV << " MeV\n";
    }
}

// ==============================================================================
// At-Rest DoIt: Cascade First, Then Capture vs Bound DIO
// ==============================================================================
void MuonCapture::atRestDoIt(Track& track, const Material& material, Step& step,
                             std::vector<Track>& secondaries) {
    const int z = static_cast<int>(std::lround(material.getZ()));
    if (z < 1) {
        track.setKineticEnergy(0.0);
        track.killTrack();
        return;
    }

    const Vector3D pos = track.getPosition();
    const double binding1s = get1sBindingEnergy(z);

    // 1. Atomic cascade always happens upon stopping
    emitAtomicCascade(z, pos, step, secondaries);

    // 2. Branch between nuclear capture and bound decay-in-orbit
    const double lambdaC = getCaptureRatePerNs(z);
    const double lambdaD = getHuffFactor(z) * getFreeDecayRatePerNs(track);
    if (RandomEngine::uniNumGen() * (lambdaC + lambdaD) < lambdaC) {
        performNuclearCapture(z, material, pos, step, secondaries);
    } else {
        performBoundDecay(binding1s, pos, step, secondaries);
    }

    track.setKineticEnergy(0.0);
    track.killTrack();
}