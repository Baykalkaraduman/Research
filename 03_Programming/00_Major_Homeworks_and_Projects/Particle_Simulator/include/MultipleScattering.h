#pragma once

#include "PhysicsProcess.h"

// ==============================================================================
// Derived Class: MultipleScattering
// Many small-angle Coulomb deflections on nuclei, added up over a step.
// Conceptually similar to G4MuMultipleScattering / G4eMultipleScattering,
// with the Highland formula instead of Geant4's Urban / WentzelVI models.
// ==============================================================================
class MultipleScattering : public PhysicsProcess {
private:
    // The log term makes the Highland formula NON-ADDITIVE: N thin steps, each with
    // theta0(x/N), added in quadrature, give LESS than theta0(x) for the whole layer
    // (e.g. 100 steps through 1 X0: ~17 % less in theta0^2). The result would
    // depend on the step size.
    //
    // Fix: remember the path length already travelled in this material, L. A step
    // of length s then adds the variance
    //       theta0^2(L + s) - theta0^2(L)
    // so the total over a layer is theta0^2 of the whole layer, for ANY steps.
    // (Reset when the material changes, or for a new track.)
    const Material* currentMaterial = nullptr;
    double pathLengthInMaterial = 0.0;

    // Highland theta0^2 for a total path `length` (-> MultipleScattering.cpp)
    // (G4UrbanMscModel::ComputeTheta0)
    double computeTheta0Squared(const Track& track, const Material& material, double length) const;

public:
    MultipleScattering() : PhysicsProcess("msc") {}

    bool isApplicable(const Particle& particle) const override {return particle.getCharge() != 0.0;}

    void startTracking(Track&) override {currentMaterial = nullptr; pathLengthInMaterial = 0.0;}

    // Does not limit the step (kInfinity, the default of PhysicsProcess). The
    // Highland formula with the path-length bookkeeping above works for any step.

    // Change the direction at the end of the step (-> MultipleScattering.cpp)
    void alongStepDoIt(Track& track, const Material& material, double stepLength) override;

    // [INFO] Public for the validation program: theta0 for one layer of thickness x
    double computeTheta0(const Track& track, const Material& material, double thickness) const;
};
