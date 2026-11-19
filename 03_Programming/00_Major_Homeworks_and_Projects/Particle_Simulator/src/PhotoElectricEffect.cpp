#include "PhotoElectricEffect.h"
#include "PhotonDataTable.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::m_e;

// ==============================================================================
// mu/rho from the XCOM table
// ==============================================================================
double PhotoElectricEffect::computeMassAttenuation(double energy, const Material& material) const {
    if (material.getZ() < 1.0) {
        return 0.0; // vacuum
    }
    return PhotonDataTable::get(static_cast<int>(std::lround(material.getZ()))).getPhotoElectric(energy);
}

double PhotoElectricEffect::getMeanFreePath(const Track& track, const Material& material) const {
    const double mu = computeMassAttenuation(track.getKineticEnergy(), material) * material.getDensity();
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// The innermost shell the photon can ionise.
double PhotoElectricEffect::getBindingEnergy(double energy, const Material& material) const {
    const PhotonDataTable& table = PhotonDataTable::get(static_cast<int>(std::lround(material.getZ())));
    double binding = 0.0;
    for (double shell : table.getShellEnergies()) {
        if (shell < energy) {
            binding = shell; // ascending: the last one below the photon energy
        }
    }
    return binding;
}

// ==============================================================================
// Sauter-Gavrila angular distribution (K-shell electrons, fully relativistic)
// ==============================================================================
Vector3D PhotoElectricEffect::samplePhotoElectronDirection(double electronKineticEnergy,
                                                           const Vector3D& photonDirection) const {
    const double energy = std::max(electronKineticEnergy, 1.0 * Units::eV);
    if (energy > 100.0 * Units::MeV) {
        return photonDirection;
    }
    const double tau = energy / m_e;
    const double gamma = 1.0 + tau;
    const double beta = std::sqrt(tau * (tau + 2.0)) / gamma;
    const double ac = (1.0 - beta) / beta;
    const double a1 = 0.5 * beta * gamma * tau * (gamma - 2.0);
    const double a2 = ac + 2.0;
    const double gtmax = 2.0 * (a1 + 1.0 / ac);

    double tsam, gtr;
    do {
        const double rand = RandomEngine::uniNumGen();
        tsam = 2.0 * ac * (2.0 * rand + a2 * std::sqrt(rand)) / (a2 * a2 - 4.0 * rand);
        gtr = (2.0 - tsam) * (a1 + 1.0 / (ac + tsam));
    } while (RandomEngine::uniNumGen() * gtmax > gtr);

    const double cosTheta = 1.0 - tsam;
    const double sinTheta = std::sqrt(std::max(0.0, tsam * (2.0 - tsam)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    Vector3D direction(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    return direction.rotateUz(photonDirection);
}

// ==============================================================================
// PostStepDoIt: photon absorbed, photoelectron emitted, B deposited
// ==============================================================================
void PhotoElectricEffect::postStepDoIt(Track& track, const Material& material, Step& step,
                                       std::vector<Track>& secondaries) {
    const double photonEnergy = track.getKineticEnergy();
    const double binding = getBindingEnergy(photonEnergy, material);
    const double electronEnergy = photonEnergy - binding;

    if (electronEnergy > 0.0) {
        const Vector3D direction = samplePhotoElectronDirection(electronEnergy, track.getDirection());
        secondaries.emplace_back(electronDef, track.getPosition(), direction, electronEnergy);
    }
    step.addTotalEnergyDeposit(binding);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] photon absorbed: e- with " << electronEnergy / Units::keV
                  << " keV, binding " << binding / Units::keV << " keV deposited\n";
    }

    track.setKineticEnergy(0.0);
    clearNumberOfInteractionLengthLeft();
    track.killTrack();
}
