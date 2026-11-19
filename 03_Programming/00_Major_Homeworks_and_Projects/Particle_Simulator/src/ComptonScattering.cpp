#include "ComptonScattering.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

// ==============================================================================
// Klein-Nishina total cross section per electron (k = E / m_e)
// ==============================================================================
// sigma = 2 pi r_e^2 { (1+k)/k^2 [ 2(1+k)/(1+2k) - ln(1+2k)/k ] + ln(1+2k)/(2k) - (1+3k)/(1+2k)^2 }
// For k -> 0 the terms cancel each other (each is ~1/k^2, the sum is the finite
// Thomson value 8/3 pi r_e^2): below k = 1e-3 (0.5 keV) the series is used.
double ComptonScattering::computeCrossSectionPerElectron(double energy) {
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

double ComptonScattering::getMeanFreePath(const Track& track, const Material& material) const {
    const double mu = material.getAtomDensity() * material.getZ()
                    * computeCrossSectionPerElectron(track.getKineticEnergy());
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// dsigma/deps = (pi r_e^2 / k) [ 1/eps + eps ] [ 1 - eps sin^2(theta) / (1 + eps^2) ]
// ==============================================================================
double ComptonScattering::computeDifferentialCrossSection(double energy, double eps) {
    const double k = energy / m_e;
    const double epsMin = 1.0 / (1.0 + 2.0 * k);
    if (eps < epsMin || eps > 1.0) {return 0.0;}
    const double oneMinusCos = (1.0 - eps) / (eps * k);
    const double sin2 = oneMinusCos * (2.0 - oneMinusCos);
    return M_PI * r_e * r_e / k * (1.0 / eps + eps) * (1.0 - eps * sin2 / (1.0 + eps * eps));
}

// ==============================================================================
// Sampling eps (Butcher & Messel 1960, as G4KleinNishinaCompton)
// ==============================================================================
// Write dsigma/deps ~ f(eps) g(eps) with
//   f(eps) = 1/eps + eps          (easy to sample: a sum of two pieces)
//   g(eps) = 1 - eps sin^2/(1 + eps^2)  in [0, 1]  (the rejection function)
// The piece 1/eps on [eps0, 1] has weight alpha1 = ln(1/eps0): eps = eps0^r.
// The piece eps has weight alpha2 = (1 - eps0^2)/2: eps^2 uniform in [eps0^2, 1].
double ComptonScattering::sampleEnergyFraction(double energy) {
    const double k = energy / m_e;
    const double eps0 = 1.0 / (1.0 + 2.0 * k);
    const double eps0sq = eps0 * eps0;
    const double alpha1 = -std::log(eps0);
    const double alpha2 = alpha1 + 0.5 * (1.0 - eps0sq);

    double eps, epsSq, greject;
    do {
        if (alpha1 > alpha2 * RandomEngine::uniNumGen()) {
            eps = std::exp(-alpha1 * RandomEngine::uniNumGen()); // eps0^r
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

// ==============================================================================
// PostStepDoIt: new photon energy and direction, recoil electron
// ==============================================================================
void ComptonScattering::postStepDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) {
    const double energy = track.getKineticEnergy();
    const double k = energy / m_e;
    const double eps = sampleEnergyFraction(energy);

    const double cosTheta = std::clamp(1.0 - (1.0 - eps) / (eps * k), -1.0, 1.0);
    const double sinTheta = std::sqrt((1.0 - cosTheta) * (1.0 + cosTheta));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    const Vector3D oldDirection = track.getDirection();
    Vector3D newDirection(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    newDirection.rotateUz(oldDirection);

    // Electron: momentum conservation, p_e = p_gamma - p_gamma'  (photon: p = E)
    const double scatteredEnergy = eps * energy;
    const double electronEnergy = energy - scatteredEnergy;
    if (electronEnergy > 0.0) {
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
