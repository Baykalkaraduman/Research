#include "Ionization.h"
#include <iostream>
#include <algorithm> // For std::max

Ionization::Ionization() : PhysicsProcess("Ionization") {}

void Ionization::doIt(Track& track, const Material& material, double stepLength, std::vector<Track>& secondaries) {
    // 1. Physical checks: Neutral particles do not suffer ionization
    if (track.getParticleDefinition()->getCharge() == 0.0) {
        return; 
    }

    // 2. Relativistic Kinematics Calculation
    double mass = track.getParticleDefinition()->getMass();
    double kinEnergy = track.getKineticEnergy();
    double betaSq = std::pow(track.getBeta(), 2);
    double gamma = track.getGamma();

    // 3. Simplified Bethe-Bloch dE/dx
    // Proportional to material density and inversely proportional to beta^2.
    // The constant 0.02 is arbitrary for our toy model to scale the GeV loss.
    double stoppingPower = (material.getDensity() * 0.02) / betaSq; 
    
    // Total energy loss in this step
    double energyLoss = stoppingPower * stepLength;

    // 4. Apply the loss to the track
    std::cout << "    [" << processName << "] Active! " 
              << "Loss: " << energyLoss << " GeV (beta^2: " << betaSq << ")\n";
              
    track.deductEnergy(energyLoss);
}