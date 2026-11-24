#pragma once

#include "PhysicalConstants.h"
#include "PhysicsProcess.h"

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
    [[nodiscard]] static double getRemainingProperTime(const Track& track);

public:
    Decay() : PhysicsProcess("Decay") {}

    // Only particles that are unstable AND have decay channels
    [[nodiscard]] bool isApplicable(const Particle& particle) const override {
        return !particle.isStable() && particle.getDecayTable() != nullptr;
    }

    [[nodiscard]] bool isAtRestApplicable(const Particle& particle, const Material&) const override {
        return isApplicable(particle);
    }

    // Draw the decay time tau * (-ln u), unless one was given already. -> Decay.cpp
    void startTracking(Track& track) override;

    // Decay length in the lab: lambda = beta gamma c tau (G4Decay::GetMeanFreePath).
    // No longer used for the stepping, kept for printing and the validation.
    [[nodiscard]] static double getMeanFreePath(const Track& track) {
        return track.getBetaGamma() * PhysicalConstants::c * track.getParticleDefinition()->getLifetime();
    }

    // In flight: distance to the decay point = remaining proper time * beta gamma c
    [[nodiscard]] double postStepGetPhysicalInteractionLength(const Track& track, const Material&) override {
        return getRemainingProperTime(track) * PhysicalConstants::c * track.getBetaGamma();
    }

    void postStepDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) override {
        decayIt(track, secondaries);
    }

    // At rest: time until the decay = remaining proper time
    [[nodiscard]] double atRestGetPhysicalInteractionLength(const Track& track, const Material&) override {
        return getRemainingProperTime(track);
    }

    void atRestDoIt(Track& track, const Material&, Step&, std::vector<Track>& secondaries) override {
        decayIt(track, secondaries);
    }
};