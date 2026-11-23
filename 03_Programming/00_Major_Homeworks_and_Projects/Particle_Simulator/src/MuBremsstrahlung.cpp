#include "MuBremsstrahlung.h"
#include "ModifiedMephi.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>

// ==============================================================================
// Constants & Helpers (Kelner, Kokoulin, Petrukhin 1995)
// ==============================================================================
namespace {

constexpr double kSqrtE = 1.6487212707001282; // sqrt(e)
constexpr double kB = 183.0;                  // nuclear screening constant
constexpr double kBPrime = 1429.0;            // atomic electron screening constant

// Maximum energy transfer on an atomic electron: k_max' = E / (1 + M^2 / (2 m_e E))
double computeKMaxPrime(double muonMass, double totalEnergy) {
    return totalEnergy / (1.0 + (muonMass * muonMass) / (2.0 * PhysicalConstants::m_e * totalEnergy));
}

} // namespace

// ==============================================================================
// EnergyLossProcess Hooks
// ==============================================================================
std::vector<double> MuBremsstrahlung::getIntegrationBreakpoints(const Particle& particle, double kineticEnergy, const Material&) const {
    const double mMu = particle.getMass();
    const double totalEnergy = kineticEnergy + mMu;
    return {computeKMaxPrime(mMu, totalEnergy)};
}

// ==============================================================================
// Differential Cross Section (dsigma / dk per atom)
// ==============================================================================
double MuBremsstrahlung::computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double transfer, const Material& material) const {
    const double z = material.getZ();
    const double a = material.getA();
    if (z < 1.0 || transfer <= 0.0 || transfer >= kineticEnergy) {
        return 0.0;
    }

    const double mMu = particle.getMass();
    const double mE = PhysicalConstants::m_e;
    const double totalEnergy = kineticEnergy + mMu;
    const double v = transfer / totalEnergy;
    const double delta = 0.5 * mMu * mMu * v / (totalEnergy - transfer);

    // Nuclear screening and finite nuclear size (Phi_n)
    const double z13 = std::cbrt(z);
    const double dN = 1.54 * std::pow(a, 0.27);
    const double dNPrime = std::pow(dN, 1.0 - 1.0 / z);
    const double bZ13 = kB / z13;

    const double numN = bZ13 * (mMu + delta * (dNPrime * kSqrtE - 2.0));
    const double denN = dNPrime * (mE + delta * kSqrtE * bZ13);
    const double phiN = (numN > denN && denN > 0.0) ? std::log(numN / denN) : 0.0;

    // Atomic electron contribution (Phi_e), active only below k_max'
    double phiE = 0.0;
    const double kMaxPrime = computeKMaxPrime(mMu, totalEnergy);
    if (transfer < kMaxPrime) {
        const double bPrimeZ23 = kBPrime / (z13 * z13);
        const double numE = bPrimeZ23 * mMu;
        const double denE = (1.0 + delta * mMu / (mE * mE * kSqrtE)) * (mE + delta * kSqrtE * bPrimeZ23);
        phiE = (numE > denE && denE > 0.0) ? std::log(numE / denE) : 0.0;
    }

    // Classical muon radius: r_mu = r_e * (m_e / M)
    const double rMu = PhysicalConstants::r_e * (mE / mMu);
    const double factor = (16.0 / 3.0) * PhysicalConstants::alpha * rMu * rMu * z;
    const double spectral = 1.0 - v + 0.75 * v * v;

    return factor * spectral * (z * phiN + phiE) / transfer;
}

// ==============================================================================
// Secondary Production
// ==============================================================================
// Samples the photon polar angle via ModifiedMephi, rotates it into the muon
// frame, emits the photon, and updates the muon energy and recoil direction.
void MuBremsstrahlung::produceSecondaries(Track& track, const Material& /*material*/, Step& /*step*/, std::vector<Track>& secondaries, double transfer) {
    if (gammaDef == nullptr) {
        return;
    }

    const Particle& muon = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();
    const double mMu = muon.getMass();

    // Sample photon direction in the local frame (z along incoming muon)
    const double cosTheta = ModifiedMephi::sampleCosTheta(kineticEnergy, transfer, mMu);
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    // Rotate photon direction relative to the incoming muon direction
    const Vector3D muDir = track.getDirection();
    Vector3D gammaDir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    gammaDir.rotateUz(muDir);
    secondaries.emplace_back(gammaDef, track.getPosition(), gammaDir, transfer);

    // Update muon energy and recoil direction (momentum conservation)
    const double pOld = track.getMomentumMag();
    const Vector3D newMomentum = muDir * pOld - gammaDir * transfer;
    track.deductEnergy(transfer);
    if (track.getKineticEnergy() > 0.0) {
        track.setDirection(newMomentum);
    }
}