#pragma once

#include "DiscreteProcess.h"
#include "PhysicalConstants.h"

// ==============================================================================
// Derived Class: Decay
// Simulates the spontaneous decay of unstable particles in flight.
// Equivalent to G4Decay (in-flight part).
// ==============================================================================
// Derives from DiscreteProcess instead of PhysicsProcess:
// a decay happens at a point, so it is a post-step process.
class Decay : public DiscreteProcess {
private:
    const Particle* decayProductDef; // The particle type this decays into (e.g., Electron)

protected:
    // Decay length in the lab:  lambda = beta * gamma * c * tau = (p / m) * c * tau
    // (beta*gamma = p/m; time dilation is inside gamma). Does not depend on the material.
    double getMeanFreePath(const Track& track, const Material&) const override {
        return track.getMomentumMag() / track.getMass() * PhysicalConstants::c * track.getParticleDefinition()->getLifetime();
    }

public:
    // One-line constructor
    explicit Decay(const Particle* productDef) : DiscreteProcess("Decay"), decayProductDef(productDef) {}

    // Tests whether the particle can decay or not
    bool isApplicable(const Particle& particle) const override {return !particle.isStable();}

    // The decay itself: produce the secondary and kill the parent
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
