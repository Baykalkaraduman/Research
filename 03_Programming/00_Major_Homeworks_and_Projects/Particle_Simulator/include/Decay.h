#pragma once

#include "PhysicsProcess.h"
#include "PhysicalConstants.h"

// ==============================================================================
// Derived Class: Decay
// Simulates the spontaneous decay of unstable particles, in flight and at rest.
// ==============================================================================
class Decay : public PhysicsProcess {
private:
    // Choose a channel, build the daughters, boost them, kill the parent.
    // Shared by the in-flight and the at-rest DoIt. -> Decay.cpp
    void decayIt(Track& track, std::vector<Track>& secondaries);

    // Proper time still left before the decay (>= 0). -> Decay.cpp
    static double getRemainingProperTime(const Track& track);

public:
    Decay() : PhysicsProcess("Decay") {}

    // Only particles that are unstable AND have decay channels
    bool isApplicable(const Particle& particle) const override { return !particle.isStable() && particle.getDecayTable() != nullptr; }

    // Draw the decay time tau * (-ln u), unless one was given already. -> Decay.cpp
    void startTracking(Track& track) override;

    // Decay length in the lab: lambda = beta gamma c tau (G4Decay::GetMeanFreePath).
    // No longer used for the stepping, kept for printing and the validation.
    static double getMeanFreePath(const Track& track) { return track.getBetaGamma() * PhysicalConstants::c * track.getParticleDefinition()->getLifetime(); }

    // In flight: distance to the decay point = remaining proper time * beta gamma c
    double postStepGetPhysicalInteractionLength(const Track& track, const Material&) override {return getRemainingProperTime(track) * PhysicalConstants::c * track.getBetaGamma();}
    void postStepDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) override {decayIt(track, secondaries);}

    // At rest: time until the decay = remaining proper time
    double atRestGetPhysicalInteractionLength(const Track& track) override {return getRemainingProperTime(track);}
    void atRestDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) override {decayIt(track, secondaries);}
};
