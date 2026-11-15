#include "Track.h"
#include <iostream>

// ==============================================================================
// Constructor Implementation
// ==============================================================================
Track::Track(const Particle* pDef, const Vector3D& initPos, const Vector3D& initDirection, double initKinEnergy)
    : particleDefinition(pDef), 
      position(initPos), 
      direction(initDirection.unit()), 
      kineticEnergy(initKinEnergy), 
      globalTime(0.0), 
      isAlive(true) 
{
    // The track starts at global time t = 0.0 with the given initial kinematics.
}

// ==============================================================================
// Setters / Kinematic Updaters
// ==============================================================================
void Track::deductEnergy(double energyLoss) {
    kineticEnergy -= energyLoss;
    // If kinetic energy drops to zero or below, the particle stops.
    if (kineticEnergy <= 0.0) {
        kineticEnergy = 0.0;
        killTrack();
    }
}

// ==============================================================================
// Display Utility
// ==============================================================================
void Track::displayState() const {
    std::cout << "  [Track State] Time: " << globalTime << " ns | Alive: " << (isAlive ? "Yes" : "No") << "\n";
    std::cout << "  -> Position : "; position.display(); std::cout << " cm\n";
    std::cout << "  -> Direction : "; direction.display(); std::cout << "\n";
    std::cout << "  -> KinEnergy: " << kineticEnergy << " GeV\n";
    std::cout << "  --------------------------------------------------\n";
}