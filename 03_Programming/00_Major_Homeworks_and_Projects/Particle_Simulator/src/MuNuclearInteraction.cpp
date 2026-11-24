#include "MuNuclearInteraction.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>

// ==============================================================================
// Constants & Helpers (Bezrukov & Bugaev 1981)
// ==============================================================================
namespace {

constexpr double kM1Sq = 0.54 * Units::GeV * Units::GeV;
constexpr double kM2Sq = 1.80 * Units::GeV * Units::GeV;

// Nuclear shadowing factor G(x) with Taylor branch for small x
double computeShadowingG(double a, double sigmaGammaNMicrobarn) {
    if (a <= 1.5) {
        return 1.0;
    }

    const double x = 0.00282 * std::cbrt(a) * sigmaGammaNMicrobarn;
    if (x < 1.0e-3) {
        return 1.0 - 0.375 * x + 0.1 * x * x;
    }

    const double x2 = x * x;
    const double x3 = x2 * x;
    return (3.0 / x3) * (0.5 * x2 - 1.0 + (1.0 + x) * std::exp(-x));
}

// Two-body CM breakup momentum for parent mass M -> m1 + m2
double twoBodyMomentum(double parentMass, double m1, double m2) {
    if (parentMass <= m1 + m2) {
        return 0.0;
    }
    const double s = parentMass * parentMass;
    const double sum = m1 + m2;
    const double diff = m1 - m2;
    return std::sqrt(std::max(0.0, (s - sum * sum) * (s - diff * diff))) / (2.0 * parentMass);
}

// Isotropic unit vector in 3D
Vector3D sampleIsotropicDirection() {
    const double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    return Vector3D(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

// Lorentz-boosts 3-momentum pStar and energy eStar from a frame moving with
// velocity betaVec and Lorentz factor gamma into the parent frame.
void boostMomentum(const Vector3D& pStar, double eStar, const Vector3D& betaVec, double gamma, Vector3D& pOut, double& eOut) {
    const double beta2 = betaVec.dot(betaVec);
    const double betaDotP = betaVec.dot(pStar);
    const double gammaFactor = (beta2 > 1.0e-16) ? ((gamma - 1.0) / beta2) * betaDotP + gamma * eStar : eStar;
    pOut = pStar + betaVec * gammaFactor;
    eOut = gamma * (eStar + betaDotP);
}

// Boosts a secondary particle from the CM frame into the lab frame and appends it.
void emitBoostedSecondary(const Particle* def, const Vector3D& position, const Vector3D& pStar, double eStar,
                          const Vector3D& betaVec, double gamma, std::vector<Track>& secondaries) {
    if (def == nullptr) {
        return;
    }

    Vector3D pLab;
    double eLab = 0.0;
    boostMomentum(pStar, eStar, betaVec, gamma, pLab, eLab);

    const double kinEnergy = std::max(0.0, eLab - def->getMass());
    const Vector3D dir = (pLab.dot(pLab) > 0.0) ? pLab.unit() : Vector3D(0.0, 0.0, 1.0);
    secondaries.emplace_back(def, position, dir, kinEnergy);
}

// Generates a charge-conserving 2-body decay (W -> p + pi) in the CM frame.
void emitTwoBodyCM(const Particle* nucleonDef, const Particle* pionDef, double w, const Vector3D& pos,
                   const Vector3D& betaCM, double gammaCM, std::vector<Track>& secondaries) {
    if (nucleonDef == nullptr || pionDef == nullptr) {
        return;
    }

    const double mN = nucleonDef->getMass();
    const double mPi = pionDef->getMass();
    const double pStarMag = twoBodyMomentum(w, mN, mPi);
    const double eNStar = std::sqrt(pStarMag * pStarMag + mN * mN);
    const double ePiStar = std::sqrt(pStarMag * pStarMag + mPi * mPi);

    const Vector3D dirStar = sampleIsotropicDirection();
    emitBoostedSecondary(nucleonDef, pos, dirStar * pStarMag, eNStar, betaCM, gammaCM, secondaries);
    emitBoostedSecondary(pionDef, pos, dirStar * (-pStarMag), ePiStar, betaCM, gammaCM, secondaries);
}

// Generates a charge-conserving 3-body phase-space decay (W -> p + pi+ + pi-) in the CM frame.
void emitThreeBodyCM(const Particle* nucleonDef, const Particle* pi1Def, const Particle* pi2Def, double w,
                     const Vector3D& pos, const Vector3D& betaCM, double gammaCM, std::vector<Track>& secondaries) {
    if (nucleonDef == nullptr || pi1Def == nullptr || pi2Def == nullptr) {
        return;
    }

    const double m1 = nucleonDef->getMass();
    const double m2 = pi1Def->getMass();
    const double m3 = pi2Def->getMass();
    const double m23Min = m2 + m3;
    const double m23Max = w - m1;
    if (m23Max <= m23Min) {
        return;
    }

    // Sample intermediate dipion mass m23 weighted by phase space p1(m23) * p2(m23)
    const double weightMax = twoBodyMomentum(w, m1, m23Min) * twoBodyMomentum(m23Max, m2, m3);
    double m23 = 0.5 * (m23Min + m23Max);
    if (weightMax > 0.0) {
        double weight = 0.0;
        do {
            m23 = m23Min + (m23Max - m23Min) * RandomEngine::uniNumGen();
            weight = twoBodyMomentum(w, m1, m23) * twoBodyMomentum(m23, m2, m3);
        } while (weight < weightMax * RandomEngine::uniNumGen());
    }

    // 1. Decay W -> nucleon (1) + dipion (23) in CM frame
    const double p1Mag = twoBodyMomentum(w, m1, m23);
    const double e1Star = std::sqrt(p1Mag * p1Mag + m1 * m1);
    const double e23Star = std::sqrt(p1Mag * p1Mag + m23 * m23);
    const Vector3D p1Star = sampleIsotropicDirection() * p1Mag;
    const Vector3D p23Star = p1Star * (-1.0);

    // 2. Decay dipion (23) -> pi1 (2) + pi2 (3) in dipion rest frame, then boost to CM
    const double p2Mag = twoBodyMomentum(m23, m2, m3);
    const double e2Dipion = std::sqrt(p2Mag * p2Mag + m2 * m2);
    const double e3Dipion = std::sqrt(p2Mag * p2Mag + m3 * m3);
    const Vector3D p2Dipion = sampleIsotropicDirection() * p2Mag;
    const Vector3D p3Dipion = p2Dipion * (-1.0);

    const Vector3D beta23 = p23Star * (1.0 / e23Star);
    const double gamma23 = e23Star / m23;
    Vector3D p2Star, p3Star;
    double e2Star = 0.0, e3Star = 0.0;
    boostMomentum(p2Dipion, e2Dipion, beta23, gamma23, p2Star, e2Star);
    boostMomentum(p3Dipion, e3Dipion, beta23, gamma23, p3Star, e3Star);

    // 3. Boost all three hadrons from CM frame to lab frame
    emitBoostedSecondary(nucleonDef, pos, p1Star, e1Star, betaCM, gammaCM, secondaries);
    emitBoostedSecondary(pi1Def, pos, p2Star, e2Star, betaCM, gammaCM, secondaries);
    emitBoostedSecondary(pi2Def, pos, p3Star, e3Star, betaCM, gammaCM, secondaries);
}

} // namespace

// ==============================================================================
// EnergyLossProcess Hooks
// ==============================================================================
double MuNuclearInteraction::getMaxEnergyTransfer(const Particle& particle, double kineticEnergy, const Material&) const {
    const double mMu = particle.getMass();
    const double mP = PhysicalConstants::m_p;
    const double totalEnergy = kineticEnergy + mMu;
    const double maxByRecoil = totalEnergy - 0.5 * mP * (1.0 + (mMu * mMu) / (mP * mP));
    return std::max(0.0, std::min(kineticEnergy, maxByRecoil));
}

// ==============================================================================
// Differential Cross Section (dsigma / deps per atom)
// ==============================================================================
double MuNuclearInteraction::computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double transfer, const Material& material) const {
    const double z = material.getZ();
    const double a = material.getA();
    if (z < 1.0 || a < 1.0 || transfer < minThreshold || transfer >= getMaxEnergyTransfer(particle, kineticEnergy, material)) {
        return 0.0;
    }

    const double mMu = particle.getMass();
    const double mMuSq = mMu * mMu;
    const double totalEnergy = kineticEnergy + mMu;
    const double v = transfer / totalEnergy;
    if (v <= 0.0 || v >= 1.0) {
        return 0.0;
    }

    const double kappa = 1.0 - 2.0 / v + 2.0 / (v * v);
    const double t = (mMuSq * v * v) / (1.0 - v);

    // Real photon-nucleon cross section sigma_gammaN(eps) in cm^2
    const double epsGeV = transfer / Units::GeV;
    const double lnTerm = std::log(0.0213 * epsGeV);
    const double sigmaGammaNUb = 114.3 + 1.647 * lnTerm * lnTerm;
    const double sigmaGammaN = sigmaGammaNUb * 1.0e-30 * Units::cm * Units::cm;

    const double shadowingG = computeShadowingG(a, sigmaGammaNUb);

    const double term1 = 0.75 * shadowingG * (kappa * std::log1p(kM1Sq / t) - (kappa * kM1Sq) / (kM1Sq + t) - (2.0 * mMuSq) / t);
    const double term2 = 0.25 * (kappa * std::log1p(kM2Sq / t) - (2.0 * mMuSq) / t);

    const double ratio2 = t / kM2Sq;
    const double logTerm2 = (ratio2 < 1.0e-6) ? (1.0 - 0.5 * ratio2) : (std::log1p(ratio2) / ratio2);
    const double term3 = (0.5 * mMuSq / t) * (0.75 * shadowingG * (kM1Sq / (kM1Sq + t)) + 0.25 * logTerm2);

    const double bracket = term1 + term2 + term3;
    if (bracket <= 0.0) {
        return 0.0;
    }

    return (PhysicalConstants::alpha / (2.0 * M_PI)) * a * sigmaGammaN * (v * v / transfer) * bracket;
}

// ==============================================================================
// Secondary Production
// ==============================================================================
// Samples virtual-photon Q^2 within both leptonic and hadronic kinematic bounds,
// scatters the muon, and decays the hadronic system (virtual photon + target
// nucleon at rest) in a charge-conserving channel before boosting to the lab frame.
void MuNuclearInteraction::produceSecondaries(Track& track, const Material& material, Step& /*step*/, std::vector<Track>& secondaries, double transfer) {
    const Particle& muon = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();
    const double mMu = muon.getMass();
    const double mP = PhysicalConstants::m_p;
    const double mPi = PhysicalConstants::m_pi;
    const double totalEnergy = kineticEnergy + mMu;

    // 1. Sample Q^2 analytically within min(hadronic, leptonic) upper limit
    const double newKineticEnergy = std::max(0.0, kineticEnergy - transfer);
    const double newTotalEnergy = newKineticEnergy + mMu;
    const double pOld = track.getMomentumMag();
    const double pNew = std::sqrt(newKineticEnergy * (newKineticEnergy + 2.0 * mMu));

    const double eePlusPpMinusM2 = std::max(totalEnergy * newTotalEnergy + pOld * pNew - mMu * mMu, 1.0e-12);
    const double q2Min = (2.0 * mMu * mMu * transfer * transfer) / eePlusPpMinusM2;
    const double q2MaxLep = 2.0 * (totalEnergy * newTotalEnergy + pOld * pNew - mMu * mMu);
    const double q2MaxHad = 2.0 * mP * (transfer - mPi) - mPi * mPi;
    const double q2Max = std::min(q2MaxHad, q2MaxLep);

    double q2 = q2Min;
    if (q2Max > q2Min) {
        const double yMin = q2Min / (q2Min + kM1Sq);
        const double yMax = q2Max / (q2Max + kM1Sq);
        const double y = yMin * std::pow(yMax / yMin, RandomEngine::uniNumGen());
        q2 = kM1Sq * y / std::max(1.0 - y, 1.0e-12);
    }

    // 2. Scattered muon kinematics
    double cosThetaMu = 1.0;
    if (pOld > 0.0 && pNew > 0.0) {
        cosThetaMu = (totalEnergy * newTotalEnergy - mMu * mMu - 0.5 * q2) / (pOld * pNew);
        cosThetaMu = std::min(1.0, std::max(-1.0, cosThetaMu));
    }
    const double sinThetaMu = std::sqrt(std::max(0.0, (1.0 - cosThetaMu) * (1.0 + cosThetaMu)));
    const double phiMu = 2.0 * M_PI * RandomEngine::uniNumGen();

    const Vector3D muDir = track.getDirection();
    Vector3D newMuDir(sinThetaMu * std::cos(phiMu), sinThetaMu * std::sin(phiMu), cosThetaMu);
    newMuDir.rotateUz(muDir);

    // Virtual photon momentum vector q = p_old - p_new
    const Vector3D qVec = muDir * pOld - newMuDir * pNew;

    track.deductEnergy(transfer);
    if (track.getKineticEnergy() > 0.0) {
        track.setDirection(newMuDir);
    }

    // 3. Charge-conserving hadronic final state in CM frame of (gamma* + target N)
    const double eHad = transfer + mP;
    const double w2 = eHad * eHad - qVec.dot(qVec);
    const double minHadMass = mP + mPi;
    if (w2 <= minHadMass * minHadMass) {
        return;
    }

    const double w = std::sqrt(w2);
    const Vector3D betaCM = qVec * (1.0 / eHad);
    const double gammaCM = eHad / w;
    const Vector3D pos = track.getPosition();

    // Select target proton (prob Z/A, Q = +1) vs target neutron (prob (A-Z)/A, Q = 0)
    const bool targetIsProton = (RandomEngine::uniNumGen() * material.getA()) < material.getZ();

    if (!targetIsProton && protonDef != nullptr && piMinusDef != nullptr) {
        // gamma* + n -> p + pi- (charge: 0 -> +1 - 1 = 0)
        emitTwoBodyCM(protonDef, piMinusDef, w, pos, betaCM, gammaCM, secondaries);
    } else if (w > mP + 2.0 * mPi && protonDef != nullptr && piPlusDef != nullptr && piMinusDef != nullptr) {
        // gamma* + p -> p + pi+ + pi- (charge: +1 -> +1 + 1 - 1 = +1)
        emitThreeBodyCM(protonDef, piPlusDef, piMinusDef, w, pos, betaCM, gammaCM, secondaries);
    } else if (protonDef != nullptr && pi0Def != nullptr) {
        // gamma* + p -> p + pi0 (charge: +1 -> +1 + 0 = +1)
        emitTwoBodyCM(protonDef, pi0Def, w, pos, betaCM, gammaCM, secondaries);
    }
}