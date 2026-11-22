#pragma once

#include "PhysicsProcess.h"

// ==============================================================================
// Abstract Class: DiscreteProcess
// Base for every process that happens at a point: decay now, and later Compton,
// photoelectric, pair production, bremsstrahlung, annihilation.
// Equivalent to G4VDiscreteProcess.
// ==============================================================================
class DiscreteProcess : public PhysicsProcess {
protected:
    // Mean free paths left before the next interaction. < 0 means "not sampled yet"
    double numberOfInteractionLengthLeft = -1.0;

    // lambda at the pre-step point of the current step
    // (G4VProcess::currentInteractionLength)
    double currentInteractionLength = GeoDefs::kInfinity;

    // Mean free path of this process for this track in this material (in cm).
    // Return kInfinity if the process cannot happen here.
    virtual double getMeanFreePath(const Track& track, const Material& material) const = 0;

    // Forget the sampled value: a new one is drawn at the next GPIL. Called after
    // the interaction happened (G4VProcess::ClearNumberOfInteractionLengthLeft)
    void clearNumberOfInteractionLengthLeft() {numberOfInteractionLengthLeft = -1.0;}

public:
    explicit DiscreteProcess(std::string name) : PhysicsProcess(std::move(name)) {}

    // A new track must not inherit the previous track's sampled value
    void startTracking(Track&) override {clearNumberOfInteractionLengthLeft();}

    // Longer than one line -> DiscreteProcess.cpp
    double postStepGetPhysicalInteractionLength(const Track& track, const Material& material) override;
    void subtractNumberOfInteractionLengthLeft(double stepLength) override;
};
