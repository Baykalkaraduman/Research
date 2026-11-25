#include "GammaConversion.h"
#include "LPMFunctions.h"
#include "ModifiedTsai.h"
#include "PhotonDataTable.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::alpha;
using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

// ==============================================================================
// Pieces of the Bethe-Heitler Formula (+ Optional LPM Factor)
// ==============================================================================
double GammaConversion::computeF(double energy, const Material& material) {
    double f = 8.0 / 3.0 * std::log(material.getZ());
    if (energy > 50.0 * Units::MeV) {
        f += 8.0 * material.computeCoulombCorrection();
    }
    return f;
}

double GammaConversion::computeDifferentialCrossSection(double energy, double eps, const Material& material,
                                                        bool enableLPM) {
    const double eps0 = m_e / energy;
    if (eps <= eps0 || eps >= 1.0 - eps0 || material.getZ() < 1.0) {
        return 0.0;
    }
    const double z = material.getZ();
    const double f = computeF(energy, material);
    const double delta = 136.0 * eps0 / std::cbrt(z) / (eps * (1.0 - eps));

    const double f1 = screenFunction1(delta);
    const double f2 = screenFunction2(delta);
    const double phi1 = (f1 + 2.0 * f2) / 6.0;
    const double phi2 = (3.0 * f2 - 1.5 * f1) / 3.0; // = Phi1 for delta > 1

    const double lRad = std::log(184.15 / std::cbrt(z));
    const double lRadPrime = std::log(1194.0 / std::cbrt(z * z));
    const double xi = lRadPrime / (lRad - material.computeCoulombCorrection());

    const double u = eps * (1.0 - eps);
    const double value = (1.0 - 2.0 * u) * (phi1 - 0.5 * f) + (2.0 / 3.0) * u * (phi2 - 0.5 * f);
    const double bh = alpha * r_e * r_e * z * (z + xi) * std::max(value, 0.0);

    if (!enableLPM || bh <= 0.0) {
        return bh;
    }
    return bh * LPMFunctions::computePairProductionFactor(energy, eps, material);
}

// sigma(E) = 2 * integral_{eps0}^{1/2} (dsigma/deps) deps.
// Splits at epsSplit = 0.05: log-spaced Simpson below (resolving narrow PeV LPM peaks
// near eps ~ E_LPM / (8 E)), linear Simpson above.
double GammaConversion::computeCrossSectionPerAtom(double energy, const Material& material, bool enableLPM) {
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return 0.0;
    }

    const double eps0 = m_e / energy;
    const double epsLow = eps0 * (1.0 + 1.0e-12);
    if (epsLow >= 0.5) {
        return 0.0;
    }

    constexpr double kEpsSplit = 0.05;
    double halfIntegral = 0.0;

    // 1. Log-spaced Simpson on [epsLow, kEpsSplit] when epsLow < kEpsSplit
    if (epsLow < kEpsSplit) {
        constexpr int nLog = 128;
        const double lnMin = std::log(epsLow);
        const double h = (std::log(kEpsSplit) - lnMin) / nLog;
        double sum = 0.0;
        for (int i = 0; i <= nLog; ++i) {
            const double w = (i == 0 || i == nLog) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
            const double eps = std::exp(lnMin + i * h);
            sum += w * eps * computeDifferentialCrossSection(energy, eps, material, enableLPM);
        }
        halfIntegral += sum * h / 3.0;
    }

    // 2. Linear Simpson on [max(epsLow, kEpsSplit), 0.5]
    const double linMin = std::max(epsLow, kEpsSplit);
    constexpr int nLin = 128;
    const double h = (0.5 - linMin) / nLin;
    double sum = 0.0;
    for (int i = 0; i <= nLin; ++i) {
        const double w = (i == 0 || i == nLin) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        sum += w * computeDifferentialCrossSection(energy, linMin + i * h, material, enableLPM);
    }
    halfIntegral += sum * h / 3.0;

    return 2.0 * halfIntegral;
}

// ==============================================================================
// Cross Section & Mean Free Path (XCOM + Optional LPM Scaling Table)
// ==============================================================================
double GammaConversion::getCrossSectionPerAtom(double energy, const Material& material, bool enableLPM) {
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return 0.0;
    }
    const double muOverRho = PhotonDataTable::get(static_cast<int>(std::lround(material.getZ()))).getPair(energy);
    const double sigmaXcom = muOverRho * material.getA() / PhysicalConstants::N_A;
    if (!enableLPM) {
        return sigmaXcom;
    }

    // Below E_LPM / 50, s'(0.5) > 5 and LPM suppression is identically 1.0
    const double eLpm = LPMFunctions::computeLPMEnergy(material);
    if (eLpm <= 0.0 || energy < eLpm / 50.0) {
        return sigmaXcom;
    }

    const double sigmaBH = computeCrossSectionPerAtom(energy, material, false);
    if (sigmaBH <= 0.0) {
        return sigmaXcom;
    }
    const double sigmaLpm = computeCrossSectionPerAtom(energy, material, true);
    return sigmaXcom * (sigmaLpm / sigmaBH);
}

const GammaConversion::LpmPairTable& GammaConversion::getLpmTable(const Material& material) const {
    auto found = lpmTables.find(&material);
    if (found != lpmTables.end()) {
        return found->second;
    }

    LpmPairTable table;
    const double eLpm = LPMFunctions::computeLPMEnergy(material);
    table.energyMin = std::max(2.0 * m_e, eLpm / 50.0);
    const int nBins = std::max(10, static_cast<int>(std::lround(
        binsPerDecade * std::log10(tableEnergyMax / table.energyMin))));
    table.lnEnergyMin = std::log(table.energyMin);
    table.dLnEnergy = std::log(tableEnergyMax / table.energyMin) / nBins;
    table.sigma.resize(static_cast<std::size_t>(nBins + 1));

    for (int i = 0; i <= nBins; ++i) {
        const double e = std::exp(table.lnEnergyMin + i * table.dLnEnergy);
        table.sigma[static_cast<std::size_t>(i)] = getCrossSectionPerAtom(e, material, true);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Built LPM pair table for " << material.getName()
                  << " (E_LPM = " << eLpm / Units::TeV << " TeV)\n";
    }
    return lpmTables.emplace(&material, std::move(table)).first->second;
}

double GammaConversion::interpolateLpmTable(const LpmPairTable& table, double energy) {
    if (energy <= table.energyMin) {
        return table.sigma.front();
    }
    const double maxPos = static_cast<double>(table.sigma.size() - 1) - 1.0e-9;
    const double position = std::clamp((std::log(energy) - table.lnEnergyMin) / table.dLnEnergy, 0.0, maxPos);
    const auto bin = static_cast<std::size_t>(position);
    const double fraction = position - static_cast<double>(bin);
    return table.sigma[bin] + fraction * (table.sigma[bin + 1] - table.sigma[bin]);
}

double GammaConversion::getMeanFreePath(const Track& track, const Material& material) const {
    const double energy = track.getKineticEnergy();
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return GeoDefs::kInfinity;
    }

    double sigma = 0.0;
    if (lpmEnabled) {
        const double eLpm = LPMFunctions::computeLPMEnergy(material);
        if (eLpm > 0.0 && energy >= eLpm / 50.0) {
            sigma = (energy <= tableEnergyMax) ? interpolateLpmTable(getLpmTable(material), energy)
                                               : getCrossSectionPerAtom(energy, material, true);
        } else {
            sigma = getCrossSectionPerAtom(energy, material, false);
        }
    } else {
        sigma = getCrossSectionPerAtom(energy, material, false);
    }

    const double mu = material.getAtomDensity() * sigma;
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// Sampling eps (G4BetheHeitlerModel + Optional LPM Rejection Factor <= 1.0)
// ==============================================================================
// Returns 0.0 if below threshold or no fraction is accepted within kMaxAttempts.
double GammaConversion::sampleEnergyFraction(double energy, const Material& material, bool enableLPM) {
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return 0.0;
    }
    const double eps0 = m_e / energy;
    if (energy < 2.0 * Units::MeV) {
        return eps0 + (0.5 - eps0) * RandomEngine::uniNumGen();
    }
    const double z = material.getZ();
    const double fz = computeF(energy, material);
    const double screenFactor = 136.0 * eps0 / std::cbrt(z);
    const double screenMax = std::exp((42.24 - fz) / 8.368) - 0.952;
    if (screenMax <= 0.0) {
        return 0.0;
    }
    const double screenMin = std::min(4.0 * screenFactor, screenMax);

    const double eps1 = 0.5 - 0.5 * std::sqrt(std::max(0.0, 1.0 - screenMin / screenMax));
    const double epsMin = std::max(eps0, eps1);
    const double epsRange = 0.5 - epsMin;
    if (epsRange <= 0.0) {
        return 0.0;
    }

    const double f10 = screenFunction1(screenMin) - fz;
    const double f20 = screenFunction2(screenMin) - fz;
    const double normF1 = std::max(f10 * epsRange * epsRange, 0.0);
    const double normF2 = std::max(1.5 * f20, 0.0);
    const double normSum = normF1 + normF2;
    if (normSum <= 0.0) {
        return 0.0;
    }

    const double eLpm = enableLPM ? LPMFunctions::computeLPMEnergy(material) : 0.0;
    const bool applyLpm = enableLPM && (eLpm > 0.0) && (energy >= eLpm / 50.0);
    const double branch1Prob = normF1 / normSum;

    constexpr int kMaxAttempts = 1000;
    for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
        double eps = 0.5;
        double greject = 0.0;
        if (branch1Prob > RandomEngine::uniNumGen()) {
            eps = 0.5 - epsRange * std::cbrt(RandomEngine::uniNumGen());
            greject = (screenFunction1(screenFactor / (eps * (1.0 - eps))) - fz) / f10;
        } else {
            eps = epsMin + epsRange * RandomEngine::uniNumGen();
            greject = (screenFunction2(screenFactor / (eps * (1.0 - eps))) - fz) / f20;
        }
        if (applyLpm) {
            greject *= LPMFunctions::computePairProductionFactor(energy, eps, material);
        }
        if (greject > RandomEngine::uniNumGen()) {
            return eps;
        }
    }
    return 0.0;
}

// ==============================================================================
// PostStepDoIt: the photon disappears, e- and e+ appear
// ==============================================================================
void GammaConversion::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    const double energy = track.getKineticEnergy();
    clearNumberOfInteractionLengthLeft();

    if (energy <= 2.0 * m_e * (1.0 + 1.0e-9) || material.getZ() < 1.0 ||
        electronDef == nullptr || positronDef == nullptr) {
        return;
    }

    const double eLpm = lpmEnabled ? LPMFunctions::computeLPMEnergy(material) : 0.0;
    const double sigma = (lpmEnabled && eLpm > 0.0 && energy >= eLpm / 50.0)
                       ? ((energy <= tableEnergyMax) ? interpolateLpmTable(getLpmTable(material), energy)
                                                     : getCrossSectionPerAtom(energy, material, true))
                       : getCrossSectionPerAtom(energy, material, false);
    if (sigma <= 0.0) {
        return;
    }

    const double eps = sampleEnergyFraction(energy, material, lpmEnabled);
    if (eps <= 0.0 || eps >= 1.0) {
        return;
    }

    double electronTotal = 0.0;
    double positronTotal = 0.0;
    if (RandomEngine::uniNumGen() > 0.5) {
        electronTotal = (1.0 - eps) * energy;
        positronTotal = eps * energy;
    } else {
        electronTotal = eps * energy;
        positronTotal = (1.0 - eps) * energy;
    }
    const double electronKinetic = std::max(electronTotal - m_e, 0.0);
    const double positronKinetic = std::max(positronTotal - m_e, 0.0);

    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    const double cosElectron = ModifiedTsai::sampleCosTheta(electronKinetic);
    const double cosPositron = ModifiedTsai::sampleCosTheta(positronKinetic);
    const double sinElectron = std::sqrt(std::max(0.0, (1.0 - cosElectron) * (1.0 + cosElectron)));
    const double sinPositron = std::sqrt(std::max(0.0, (1.0 - cosPositron) * (1.0 + cosPositron)));

    Vector3D electronDirection(sinElectron * std::cos(phi), sinElectron * std::sin(phi), cosElectron);
    Vector3D positronDirection(-sinPositron * std::cos(phi), -sinPositron * std::sin(phi), cosPositron);
    electronDirection.rotateUz(track.getDirection());
    positronDirection.rotateUz(track.getDirection());

    secondaries.emplace_back(electronDef, track.getPosition(), electronDirection, electronKinetic);
    secondaries.emplace_back(positronDef, track.getPosition(), positronDirection, positronKinetic);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] gamma " << energy / Units::MeV << " MeV -> e- "
                  << electronKinetic / Units::MeV << " MeV + e+ " << positronKinetic / Units::MeV << " MeV\n";
    }

    track.setKineticEnergy(0.0);
    track.killTrack();
}