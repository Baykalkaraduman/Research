#include "Annihilation.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

// ==============================================================================
// Heitler's formula, per electron (G4eeToTwoGammaModel::ComputeCrossSectionPerElectron)
// ==============================================================================
//   sigma = pi r_e^2 / (gamma + 1) [ (gamma^2 + 4 gamma + 1)/(gamma^2 - 1) ln(gamma + sqrt(gamma^2 - 1)) - (gamma + 3) / sqrt(gamma^2 - 1) ]
// gamma: Lorentz factor of the positron (the electron is at rest).
double Annihilation::computeCrossSectionPerElectron(double kineticEnergy) {
    const double tau = std::max(kineticEnergy, 1.0 * Units::eV) / m_e;
    const double gamma = tau + 1.0;
    const double gamma2 = gamma * gamma;
    const double root = std::sqrt(gamma2 - 1.0);
    return M_PI * r_e * r_e / (gamma + 1.0)
         * ((gamma2 + 4.0 * gamma + 1.0) / (gamma2 - 1.0) * std::log(gamma + root) - (gamma + 3.0) / root);
}

double Annihilation::getMeanFreePath(const Track& track, const Material& material) const {
    if (track.getKineticEnergy() <= 0.0) {
        return GeoDefs::kInfinity;
    }
    const double mu = material.getAtomDensity() * material.getZ() * computeCrossSectionPerElectron(track.getKineticEnergy());
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// In flight (G4eeToTwoGammaModel::SampleSecondaries)
// ==============================================================================
// eps = E1 / (T + 2 m_e), the energy share of the first photon, lies between
// 1/2 -/+ sqrt(tau/(tau+2))/2. Sampled from 1/eps (eps = epsMin (epsMax/epsMin)^r)
// and accepted with the rest of Heitler's differential cross section. The angle
// of photon 1 then follows from the kinematics; photon 2 from momentum conservation.
void Annihilation::postStepDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) {
    const double kineticEnergy = track.getKineticEnergy();
    const double tau = kineticEnergy / m_e;
    const double gamma = tau + 1.0;
    const double tau2 = tau + 2.0;
    const double sqgrate = 0.5 * std::sqrt(tau / tau2);
    const double sqg2m1 = std::sqrt(tau * tau2);

    const double epsMin = 0.5 - sqgrate;
    const double epsMax = 0.5 + sqgrate;
    const double epsRatio = epsMax / epsMin;

    double eps, greject;
    do {
        eps = epsMin * std::exp(std::log(epsRatio) * RandomEngine::uniNumGen());
        greject = 1.0 - eps + (2.0 * gamma * eps - 1.0) / (eps * tau2 * tau2);
    } while (greject < RandomEngine::uniNumGen());

    const double cosTheta = std::clamp((eps * tau2 - 1.0) / (eps * sqg2m1), -1.0, 1.0);
    const double sinTheta = std::sqrt((1.0 + cosTheta) * (1.0 - cosTheta));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    const double totalEnergy = kineticEnergy + 2.0 * m_e;
    const double energy1 = eps * totalEnergy;
    const double energy2 = totalEnergy - energy1;

    Vector3D direction1(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    direction1.rotateUz(track.getDirection());
    const Vector3D direction2 = (track.getMomentum() - direction1 * energy1).unit();

    secondaries.emplace_back(gammaDef, track.getPosition(), direction1, energy1);
    secondaries.emplace_back(gammaDef, track.getPosition(), direction2, energy2);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] in flight: e+ " << kineticEnergy / Units::MeV << " MeV -> gammas "
                  << energy1 / Units::MeV << " + " << energy2 / Units::MeV << " MeV\n";
    }
    track.setKineticEnergy(0.0);
    clearNumberOfInteractionLengthLeft();
    track.killTrack();
}

// ==============================================================================
// At rest: two 511 keV photons, back to back, isotropic
// ==============================================================================
void Annihilation::atRestDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) {
    const double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinTheta = std::sqrt((1.0 - cosTheta) * (1.0 + cosTheta));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    const Vector3D direction(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);

    secondaries.emplace_back(gammaDef, track.getPosition(), direction, m_e);
    secondaries.emplace_back(gammaDef, track.getPosition(), direction * -1.0, m_e);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] at rest: two 511 keV photons\n";
    }
    track.killTrack();
}
