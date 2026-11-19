#include "GammaConversion.h"
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
// Pieces of the Bethe-Heitler formula
// ==============================================================================
// F = (8/3) ln Z is the Thomas-Fermi screening term. Above 50 MeV Geant4 adds
// the Coulomb correction 8 f_c (the e-/e+ are fast enough for the Born
// approximation's error to matter); below, it is left out. We do the same, so
// sampling and cross section stay consistent with G4BetheHeitlerModel.
double GammaConversion::computeF(double energy, const Material& material) {
    double f = 8.0 / 3.0 * std::log(material.getZ());
    if (energy > 50.0 * Units::MeV) {
        f += 8.0 * material.computeCoulombCorrection();
    }
    return f;
}

// dsigma/deps per atom (cm^2). Screening functions with F subtracted can go
// negative far from the region where they are fitted (near threshold): the
// cross section is cut at 0 there.
double GammaConversion::computeDifferentialCrossSection(double energy, double eps, const Material& material) {
    const double eps0 = m_e / energy;
    if (eps <= eps0 || eps >= 1.0 - eps0) {
        return 0.0;
    }
    const double z = material.getZ();
    const double f = computeF(energy, material);
    const double delta = 136.0 * eps0 / std::cbrt(z) / (eps * (1.0 - eps));

    // Phi1, Phi2 back from F1 = 3 Phi1 - Phi2 and F2 = 3/2 Phi1 + 1/2 Phi2
    const double f1 = screenFunction1(delta), f2 = screenFunction2(delta);
    const double phi1 = (f1 + 2.0 * f2) / 6.0;
    const double phi2 = (3.0 * f2 - 1.5 * f1) / 3.0;  // = Phi1 for delta > 1

    // xi(Z): pair production on the Z atomic electrons (Tsai: L'rad / (Lrad - f_c))
    const double lRad = std::log(184.15 / std::cbrt(z));
    const double lRadPrime = std::log(1194.0 / std::cbrt(z * z));
    const double xi = lRadPrime / (lRad - material.computeCoulombCorrection());

    const double u = eps * (1.0 - eps);
    const double value = (1.0 - 2.0 * u) * (phi1 - 0.5 * f) + 2.0 / 3.0 * u * (phi2 - 0.5 * f);
    return alpha * r_e * r_e * z * (z + xi) * std::max(value, 0.0);
}

// sigma(E) = integral over eps. The formula is symmetric in eps <-> 1 - eps, so
// twice the integral from eps0 to 1/2 (Simpson's rule, 200 intervals).
double GammaConversion::computeCrossSectionPerAtom(double energy, const Material& material) {
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return 0.0;
    }
    const double eps0 = m_e / energy;
    const int n = 200;
    const double h = (0.5 - eps0) / n;
    double sum = computeDifferentialCrossSection(energy, eps0 + 1e-12, material)
               + computeDifferentialCrossSection(energy, 0.5, material);
    for (int i = 1; i < n; ++i) {
        sum += (i % 2 == 1 ? 4.0 : 2.0) * computeDifferentialCrossSection(energy, eps0 + i * h, material);
    }
    return 2.0 * sum * h / 3.0;
}

// ==============================================================================
// XCOM tracking used for the cross section
// ==============================================================================
double GammaConversion::getCrossSectionPerAtom(double energy, const Material& material) {
    if (energy <= 2.0 * m_e || material.getZ() < 1.0) {
        return 0.0;
    }
    const double muOverRho = PhotonDataTable::get(static_cast<int>(std::lround(material.getZ()))).getPair(energy);
    return muOverRho * material.getA() / PhysicalConstants::N_A; // cm^2/g -> cm^2 per atom
}

double GammaConversion::getMeanFreePath(const Track& track, const Material& material) const {
    const double mu = material.getAtomDensity() * getCrossSectionPerAtom(track.getKineticEnergy(), material);
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// Sampling eps (G4BetheHeitlerModel::SampleSecondaries)
// ==============================================================================
// By symmetry only eps in [epsMin, 1/2] is sampled; who gets eps is decided later.
// dsigma/deps = (eps - 1/2)^2 F1(delta) + F2(delta)/2 (up to a constant), so:
//   piece 1: (eps - 1/2)^2 on [epsMin, 1/2] -> eps = 1/2 - range * r^(1/3), accept with F1(delta)/F1(deltaMin)
//   piece 2: flat                         -> eps = epsMin + range * r,     accept with F2(delta)/F2(deltaMin)
// (F decreases with delta, and delta is smallest at eps = 1/2: that is the maximum.)
// Below 2 MeV the shape hardly matters: eps uniform.
double GammaConversion::sampleEnergyFraction(double energy, const Material& material) {
    const double eps0 = m_e / energy;
    if (energy < 2.0 * Units::MeV) {
        return eps0 + (0.5 - eps0) * RandomEngine::uniNumGen();
    }
    const double z = material.getZ();
    const double fz = computeF(energy, material);
    const double screenFactor = 136.0 * eps0 / std::cbrt(z);
    const double screenMax = std::exp((42.24 - fz) / 8.368) - 0.952; // where F1 - FZ reaches 0
    const double screenMin = std::min(4.0 * screenFactor, screenMax);   // delta at eps = 1/2

    const double eps1 = 0.5 - 0.5 * std::sqrt(1.0 - screenMin / screenMax);
    const double epsMin = std::max(eps0, eps1);
    const double epsRange = 0.5 - epsMin;

    const double f10 = screenFunction1(screenMin) - fz;
    const double f20 = screenFunction2(screenMin) - fz;
    const double normF1 = std::max(f10 * epsRange * epsRange, 0.0);
    const double normF2 = std::max(1.5 * f20, 0.0);

    double eps, greject;
    do {
        if (normF1 / (normF1 + normF2) > RandomEngine::uniNumGen()) {
            eps = 0.5 - epsRange * std::cbrt(RandomEngine::uniNumGen());
            greject = (screenFunction1(screenFactor / (eps * (1.0 - eps))) - fz) / f10;
        } else {
            eps = epsMin + epsRange * RandomEngine::uniNumGen();
            greject = (screenFunction2(screenFactor / (eps * (1.0 - eps))) - fz) / f20;
        }
    } while (greject < RandomEngine::uniNumGen());
    return eps;
}

// ==============================================================================
// PostStepDoIt: the photon disappears, e- and e+ appear
// ==============================================================================
void GammaConversion::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    const double energy = track.getKineticEnergy();
    const double eps = sampleEnergyFraction(energy, material);

    // Which of the two gets eps: 50/50 (the formula is symmetric)
    double electronTotal, positronTotal;
    if (RandomEngine::uniNumGen() > 0.5) {
        electronTotal = (1.0 - eps) * energy;
        positronTotal = eps * energy;
    } else {
        electronTotal = eps * energy;
        positronTotal = (1.0 - eps) * energy;
    }
    const double electronKinetic = std::max(electronTotal - m_e, 0.0);
    const double positronKinetic = std::max(positronTotal - m_e, 0.0);

    // Angles: Tsai for each, opposite azimuths (transverse momenta roughly balance)
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    const double cosElectron = ModifiedTsai::sampleCosTheta(electronKinetic);
    const double cosPositron = ModifiedTsai::sampleCosTheta(positronKinetic);
    const double sinElectron = std::sqrt((1.0 - cosElectron) * (1.0 + cosElectron));
    const double sinPositron = std::sqrt((1.0 - cosPositron) * (1.0 + cosPositron));

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
    clearNumberOfInteractionLengthLeft();
    track.killTrack();
}
