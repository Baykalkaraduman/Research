#include "ComptonScattering.h"
#include "PhotoElectricShellData.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

namespace {

constexpr double kHcEvAngstrom = 12398.419843320026; // h*c [eV * Angstrom]

} // namespace

// ==============================================================================
// Cross Sections & Mean Free Path
// ==============================================================================
double ComptonScattering::computeCrossSectionPerElectron(double energy) {
    if (energy <= 0.0) {
        return 0.0;
    }
    const double k = energy / m_e;
    const double thomson = 8.0 / 3.0 * M_PI * r_e * r_e;
    if (k < 1.0e-3) {
        return thomson * (1.0 - 2.0 * k + 5.2 * k * k);
    }
    const double l = std::log(1.0 + 2.0 * k);
    const double a = 1.0 + 2.0 * k;
    return 2.0 * M_PI * r_e * r_e * ((1.0 + k) / (k * k) * (2.0 * (1.0 + k) / a - l / k)
                                     + l / (2.0 * k) - (1.0 + 3.0 * k) / (a * a));
}

double ComptonScattering::computeCrossSectionPerAtom(double energy, const Material& material) const {
    if (energy <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    const int z = static_cast<int>(std::lround(material.getZ()));
    const double barn = 1.0e-24 * Units::cm * Units::cm;
    return PhotoElectricShellData::get(z).getIncoherentCrossSection(energy) * barn;
}

double ComptonScattering::getMeanFreePath(const Track& track, const Material& material) const {
    const double mu = material.getAtomDensity() * computeCrossSectionPerAtom(track.getKineticEnergy(), material);
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// dsigma/deps = (pi r_e^2 / k) [ 1/eps + eps ] [ 1 - eps sin^2(theta) / (1 + eps^2) ]
// ==============================================================================
double ComptonScattering::computeDifferentialCrossSection(double energy, double eps) {
    const double k = energy / m_e;
    const double epsMin = 1.0 / (1.0 + 2.0 * k);
    if (eps < epsMin || eps > 1.0) {
        return 0.0;
    }
    const double oneMinusCos = (1.0 - eps) / (eps * k);
    const double sin2 = oneMinusCos * (2.0 - oneMinusCos);
    return M_PI * r_e * r_e / k * (1.0 / eps + eps) * (1.0 - eps * sin2 / (1.0 + eps * eps));
}

// ==============================================================================
// Sampling eps: Free Klein-Nishina & Bound S(x)/Z (Butcher & Messel 1960)
// ==============================================================================
double ComptonScattering::sampleEnergyFraction(double energy) {
    const double k = energy / m_e;
    const double eps0 = 1.0 / (1.0 + 2.0 * k);
    const double eps0sq = eps0 * eps0;
    const double alpha1 = -std::log(eps0);
    const double alpha2 = alpha1 + 0.5 * (1.0 - eps0sq);

    double eps = 1.0;
    double epsSq = 1.0;
    double greject = 0.0;
    do {
        if (alpha1 > alpha2 * RandomEngine::uniNumGen()) {
            eps = std::exp(-alpha1 * RandomEngine::uniNumGen());
            epsSq = eps * eps;
        } else {
            epsSq = eps0sq + (1.0 - eps0sq) * RandomEngine::uniNumGen();
            eps = std::sqrt(epsSq);
        }
        const double oneMinusCos = (1.0 - eps) / (eps * k);
        const double sin2 = oneMinusCos * (2.0 - oneMinusCos);
        greject = 1.0 - eps * sin2 / (1.0 + epsSq);
    } while (greject < RandomEngine::uniNumGen());

    return eps;
}

double ComptonScattering::sampleEnergyFraction(double energy, int z) {
    if (z < 1) {
        return sampleEnergyFraction(energy);
    }

    const PhotoElectricShellData& shellData = PhotoElectricShellData::get(z);
    const double invHc = 1.0 / (kHcEvAngstrom * Units::eV);
    const double invZ = 1.0 / static_cast<double>(z);

    // Prevent an infinite loop if S(x) vanishes across the entire accessible x range
    if (shellData.getScatteringFunction(energy * invHc) <= 0.0) {
        return 1.0;
    }

    const double k = energy / m_e;
    const double eps0 = 1.0 / (1.0 + 2.0 * k);
    const double eps0sq = eps0 * eps0;
    const double alpha1 = -std::log(eps0);
    const double alpha2 = alpha1 + 0.5 * (1.0 - eps0sq);

    double eps = 1.0;
    double epsSq = 1.0;
    double greject = 0.0;
    do {
        if (alpha1 > alpha2 * RandomEngine::uniNumGen()) {
            eps = std::exp(-alpha1 * RandomEngine::uniNumGen());
            epsSq = eps * eps;
        } else {
            epsSq = eps0sq + (1.0 - eps0sq) * RandomEngine::uniNumGen();
            eps = std::sqrt(epsSq);
        }
        const double oneMinusCos = (1.0 - eps) / (eps * k);
        const double sin2 = oneMinusCos * (2.0 - oneMinusCos);
        const double x = energy * invHc * std::sqrt(std::max(0.0, 0.5 * oneMinusCos));
        greject = (1.0 - eps * sin2 / (1.0 + epsSq)) * (shellData.getScatteringFunction(x) * invZ);
    } while (greject < RandomEngine::uniNumGen());

    return eps;
}

// ==============================================================================
// PostStepDoIt: Scattered Photon & Recoil Electron (T = E - E')
// ==============================================================================
void ComptonScattering::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    const double energy = track.getKineticEnergy();
    const double k = energy / m_e;
    const int z = static_cast<int>(std::lround(material.getZ()));
    const double eps = sampleEnergyFraction(energy, z);

    const double cosTheta = std::clamp(1.0 - (1.0 - eps) / (eps * k), -1.0, 1.0);
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    const Vector3D oldDirection = track.getDirection();
    Vector3D newDirection(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    newDirection.rotateUz(oldDirection);

    // Electron: momentum conservation, p_e = p_gamma - p_gamma'  (photon: p = E)
    const double scatteredEnergy = eps * energy;
    const double electronEnergy = energy - scatteredEnergy;
    if (electronEnergy > 0.0 && electronDef != nullptr) {
        const Vector3D electronDirection = (oldDirection * energy - newDirection * scatteredEnergy).unit();
        secondaries.emplace_back(electronDef, track.getPosition(), electronDirection, electronEnergy);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] E: " << energy / Units::MeV << " -> " << scatteredEnergy / Units::MeV
                  << " MeV, theta " << std::acos(cosTheta) << " rad, e- with " << electronEnergy / Units::MeV << " MeV\n";
    }

    track.setKineticEnergy(scatteredEnergy);
    track.setDirection(newDirection);
    clearNumberOfInteractionLengthLeft();
}