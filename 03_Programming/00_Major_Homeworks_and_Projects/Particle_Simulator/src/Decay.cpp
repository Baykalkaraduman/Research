#include "Decay.h"
#include "RandomEngine.h"

#include <iostream>
#include <cmath>

// ==============================================================================
// PostStepDoIt: only called when Decay limited the step, i.e. the particle has
// reached its sampled decay point. No probability check here any more: WHETHER
// it decays was already decided by the GPIL (DiscreteProcess).
// ==============================================================================
void Decay::postStepDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) {
    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] *** DECAY TRIGGERED! ***\n";
    }

    // Generate random angles for isotropic emission (Sphere surface)
    double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    double theta = std::acos(cosTheta);
    double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    // ==========================================================================
    // SECONDARY PARTICLE GENERATION
    // ==========================================================================
    // Assign remaining kinetic energy to the secondary (Simplified toy physics)
    double productKinEnergy = track.getKineticEnergy() * 0.8; // e.g., Neutrinos carry 20% away

    // Create momentum vector from spherical coordinates
    Vector3D productDirection = Vector3D::fromSpherical(1.0, theta, phi);

    // Create the secondary track and push to the stack
    Track secondaryTrack(decayProductDef, track.getPosition(), productDirection, productKinEnergy);
    secondaries.push_back(secondaryTrack);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Spawned Secondary: " << decayProductDef->getName()
                  << " at Theta: " << theta << ", Phi: " << phi << "\n";
    }

    // ==========================================================================
    // TERMINATE PARENT TRACK
    // ==========================================================================
    clearNumberOfInteractionLengthLeft(); // interaction happened
    track.killTrack();
}
