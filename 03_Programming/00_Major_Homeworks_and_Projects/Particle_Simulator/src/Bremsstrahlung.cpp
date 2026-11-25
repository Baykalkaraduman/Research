#include "Bremsstrahlung.h"
#include "LPMFunctions.h"
#include "ModifiedTsai.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::alpha;
using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

// ==============================================================================
// The differential cross section
// ==============================================================================
// Screening variables (G4eBremsstrahlungRelModel):
//   gamma = 100 m_e k / (E E' Z^(1/3))     nuclear field
//   eps   = 100 m_e k / (E E' Z^(2/3))     atomic-electron field
// gamma -> 0: complete screening (the atom's electrons hide the nucleus far away),
// gamma >> 1: no screening. Tsai's fits of the screening functions:
//   Phi1   = 16.863 - 2 ln(1 + 0.311877 g^2) + 2.4 e^(-0.9 g) + 1.6 e^(-1.5 g)
//   Phi1m2 = 2 / (3 (1 + 6.5 g + 6 g^2))                (Phi1 - Phi2)
//   Psi1   = 24.34 - 2 ln(1 + 13.111641 e^2) + 2.8 e^(-8 e) + 1.2 e^(-29.2 e)
//   Psi1m2 = 2 / (3 (1 + 40 e + 400 e^2))
// At g = e = 0: Phi1/4 = ln 184.15 and Psi1/4 = ln 1194, the radiation logarithms.
double Bremsstrahlung::computeDxSection(double kineticEnergy, double photonEnergy, const Material& material,
                                        bool enableLPM, bool includeDielectric) {
    const double totalEnergy = kineticEnergy + m_e;
    const double y = photonEnergy / totalEnergy;
    if (y <= 0.0 || photonEnergy > kineticEnergy) {
        return 0.0;
    }
    const double z = material.getZ();
    const double lnZ = std::log(z);
    const double z13 = std::cbrt(z);

    const double dum1 = 100.0 * m_e * y / (totalEnergy - photonEnergy); // 100 m_e k / (E E')
    const double g = dum1 / z13;
    const double e = dum1 / (z13 * z13);

    const double phi1 = 16.863 - 2.0 * std::log(1.0 + 0.311877 * g * g) + 2.4 * std::exp(-0.9 * g) + 1.6 * std::exp(-1.5 * g);
    const double phi1m2 = 2.0 / (3.0 * (1.0 + 6.5 * g + 6.0 * g * g));
    const double psi1 = 24.34 - 2.0 * std::log(1.0 + 13.111641 * e * e) + 2.8 * std::exp(-8.0 * e) + 1.2 * std::exp(-29.2 * e);
    const double psi1m2 = 2.0 / (3.0 * (1.0 + 40.0 * e + 400.0 * e * e));

    const double fz = lnZ / 3.0 + material.computeCoulombCorrection();
    const double main = (1.0 - y + 0.75 * y * y) * ((0.25 * phi1 - fz) + (0.25 * psi1 - 2.0 * lnZ / 3.0) / z);
    const double extra = 0.125 * (1.0 - y) * (phi1m2 + psi1m2 / z);
    const double bh = std::max(main + extra, 0.0);

    if (bh <= 0.0) {
        return 0.0;
    }
    if (enableLPM) {
        return bh * LPMFunctions::computeBremsstrahlungFactor(totalEnergy, photonEnergy, material, includeDielectric);
    }
    if (includeDielectric) {
        return bh * LPMFunctions::computeDielectricFactor(totalEnergy, photonEnergy, material);
    }
    return bh;
}

double Bremsstrahlung::computeDifferentialCrossSection(double kineticEnergy, double photonEnergy,
                                                       const Material& material, bool enableLPM,
                                                       bool includeDielectric) {
    if (photonEnergy <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    const double z = material.getZ();
    return 16.0 / 3.0 * alpha * r_e * r_e * z * z / photonEnergy
         * computeDxSection(kineticEnergy, photonEnergy, material, enableLPM, includeDielectric);
}

// sigma(k > kMin) = C Z^2 integral dxsec(k) dk/k = C Z^2 integral dxsec d(ln k)   (Simpson, 128 intervals)
double Bremsstrahlung::computeCrossSectionPerAtom(double kineticEnergy, double kMin, const Material& material,
                                                  bool enableLPM, bool includeDielectric) {
    if (kineticEnergy <= kMin || kMin <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    constexpr int n = 128;
    const double lnMin = std::log(kMin);
    const double h = (std::log(kineticEnergy) - lnMin) / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double weight = (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        sum += weight * computeDxSection(kineticEnergy, std::exp(lnMin + i * h), material,
                                         enableLPM, includeDielectric);
    }
    const double z = material.getZ();
    return 16.0 / 3.0 * alpha * r_e * r_e * z * z * sum * h / 3.0;
}

// integral_0^kMax k dsigma/dk dk = C Z^2 integral dxsec(k) dk   (Simpson, 128 intervals)
double Bremsstrahlung::computeEnergyLossPerAtom(double kineticEnergy, double kMax, const Material& material,
                                                bool enableLPM, bool includeDielectric) {
    const double upper = std::min(kMax, kineticEnergy);
    if (upper <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    constexpr int n = 128;
    const double h = upper / n;
    double sum = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double weight = (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        // k = 0 itself: y -> 0, take the limit with a tiny k
        sum += weight * computeDxSection(kineticEnergy, std::max(i * h, 1.0e-9 * upper), material,
                                         enableLPM, includeDielectric);
    }
    const double z = material.getZ();
    return 16.0 / 3.0 * alpha * r_e * r_e * z * z * sum * h / 3.0;
}

double Bremsstrahlung::computeRadiativeDEDX(double kineticEnergy, const Material& material,
                                            bool enableLPM, bool includeDielectric) {
    return material.getAtomDensity()
         * computeEnergyLossPerAtom(kineticEnergy, kineticEnergy, material, enableLPM, includeDielectric);
}

// ==============================================================================
// Tables
// ==============================================================================
const Bremsstrahlung::BremsTable& Bremsstrahlung::getTable(const Material& material) const {
    auto found = tables.find(&material);
    if (found != tables.end()) {
        return found->second;
    }
    BremsTable table;
    const int nBins = static_cast<int>(std::lround(binsPerDecade * std::log10(tableEnergyMax / tableEnergyMin)));
    table.lnEnergyMin = std::log(tableEnergyMin);
    table.dLnEnergy = std::log(tableEnergyMax / tableEnergyMin) / nBins;
    table.sigma.resize(static_cast<std::size_t>(nBins + 1));
    table.softLoss.resize(static_cast<std::size_t>(nBins + 1));
    for (int i = 0; i <= nBins; ++i) {
        const double energy = std::exp(table.lnEnergyMin + i * table.dLnEnergy);
        const auto idx = static_cast<std::size_t>(i);
        table.sigma[idx] = computeCrossSectionPerAtom(energy, productionCut, material, lpmEnabled, dielectricEnabled);
        table.softLoss[idx] = computeEnergyLossPerAtom(energy, productionCut, material, lpmEnabled, dielectricEnabled);
    }
    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Built tables for " << material.getName()
                  << " (k_c = " << productionCut / Units::keV << " keV, LPM=" << (lpmEnabled ? "on" : "off")
                  << ", DE=" << (dielectricEnabled ? "on" : "off") << ")\n";
    }
    return tables.emplace(&material, std::move(table)).first->second;
}

// Linear in ln(T). Below the table: 0 (1 keV electrons do not radiate anything that matters).
double Bremsstrahlung::interpolate(const BremsTable& table, const std::vector<double>& values, double energy) {
    if (energy <= tableEnergyMin) {
        return 0.0;
    }
    const double maxPos = static_cast<double>(values.size() - 1) - 1.0e-9;
    const double position = std::clamp((std::log(energy) - table.lnEnergyMin) / table.dLnEnergy, 0.0, maxPos);
    const auto bin = static_cast<std::size_t>(position);
    const double fraction = position - static_cast<double>(bin);
    return values[bin] + fraction * (values[bin + 1] - values[bin]);
}

double Bremsstrahlung::getMeanFreePath(const Track& track, const Material& material) const {
    if (material.getZ() < 1.0 || track.getKineticEnergy() <= productionCut) {
        return GeoDefs::kInfinity;
    }
    const BremsTable& table = getTable(material);
    const double mu = material.getAtomDensity() * interpolate(table, table.sigma, track.getKineticEnergy());
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// Along step: photons below the cut, a continuous loss
// ==============================================================================
// It is only k_c / E of the radiative loss (the k dsigma/dk spectrum is ~flat), so
// a few % at most; it uses the pre-step energy, like Ionization for short steps.
void Bremsstrahlung::alongStepDoIt(Track& track, const Material& material, double stepLength) {
    if (material.getZ() < 1.0) {
        return;
    }
    const BremsTable& table = getTable(material);
    const double lossPerCm = material.getAtomDensity() * interpolate(table, table.softLoss, track.getKineticEnergy());
    const double loss = std::min(lossPerCm * stepLength, track.getKineticEnergy());
    if (loss > 0.0) {
        track.deductEnergy(loss);
    }
}

// ==============================================================================
// Sampling k between kMin and T
// ==============================================================================
// dsigma/dk = (1/k) * dxsec(k). 1/k: ln k uniform in [ln kMin, ln T]. Then accept
// with dxsec(k) / dxsecMax. dxsec is largest at y -> 0 with complete screening
// (all prefactors <= 1 there, all screening functions at their maxima, and
// LPM/dielectric suppression ratios are guaranteed <= 1.0):
//   dxsecMax = (ln 184.15 - ln Z/3 - f_c) + (ln 1194 - 2/3 ln Z) / Z + (1 + 1/Z) / 12
// Returns 0.0 if no photon is accepted within kMaxAttempts (Geant4-style guard).
double Bremsstrahlung::samplePhotonEnergy(double kineticEnergy, double kMin, const Material& material,
                                          bool enableLPM, bool includeDielectric) {
    if (kineticEnergy <= kMin || kMin <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    const double z = material.getZ();
    const double lnZ = std::log(z);
    const double dxsecMax = (std::log(184.15) - lnZ / 3.0 - material.computeCoulombCorrection())
                          + (std::log(1194.0) - 2.0 * lnZ / 3.0) / z + (1.0 + 1.0 / z) / 12.0;
    if (dxsecMax <= 0.0) {
        return 0.0;
    }

    const double lnRatio = std::log(kineticEnergy / kMin);
    constexpr int kMaxAttempts = 1000;
    for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
        const double k = kMin * std::exp(lnRatio * RandomEngine::uniNumGen());
        const double dxsec = computeDxSection(kineticEnergy, k, material, enableLPM, includeDielectric);
        if (dxsec > dxsecMax * RandomEngine::uniNumGen()) {
            return k;
        }
    }
    return 0.0;
}

// ==============================================================================
// Post step: one photon above the cut
// ==============================================================================
void Bremsstrahlung::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    const double kineticEnergy = track.getKineticEnergy();
    clearNumberOfInteractionLengthLeft();

    // Along-step loss may have dropped T to/below the cut, or the discrete cross section may vanish
    if (kineticEnergy <= productionCut * (1.0 + 1.0e-9) ||
        computeCrossSectionPerAtom(kineticEnergy, productionCut, material, lpmEnabled, dielectricEnabled) <= 0.0) {
        return;
    }

    const double k = samplePhotonEnergy(kineticEnergy, productionCut, material, lpmEnabled, dielectricEnabled);
    if (k <= 0.0 || k > kineticEnergy) {
        return;
    }

    // Photon inside the 1/gamma cone of the electron
    const double cosTheta = ModifiedTsai::sampleCosTheta(kineticEnergy);
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    Vector3D photonDirection(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    photonDirection.rotateUz(track.getDirection());
    secondaries.emplace_back(gammaDef, track.getPosition(), photonDirection, k);

    // Electron: momentum conservation (the nucleus takes the rest, its recoil energy is negligible)
    const double momentum = track.getMomentumMag();
    const Vector3D newDirection = (track.getDirection() * momentum - photonDirection * k).unit();
    track.setDirection(newDirection);
    track.deductEnergy(k);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] photon of " << k / Units::MeV << " MeV, e "
                  << kineticEnergy / Units::MeV << " -> " << track.getKineticEnergy() / Units::MeV << " MeV\n";
    }
}