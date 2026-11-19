#pragma once

#include "DiscreteProcess.h"

// ==============================================================================
// Derived Class: Annihilation
// New. e+ e- -> gamma gamma with an electron of the material.
//
// Two ways, both here:
// in flight (post-step, discrete): Heitler cross section. It falls with energy
// (~ln(2 gamma)/gamma), so most positrons first slow down and stop and the stopped positron 
// annihilates immediately, two photons of m_e c^2 = 511 keV, back to back, isotropic.
// The energy released is T + 2 m_e: the positron's rest mass AND that of the
// atomic electron. That is why a shower does not lose the 2 m_e each pair production put into mass.
// ==============================================================================
class Annihilation : public DiscreteProcess {
private:
    const Particle* gammaDef;

protected:
    // lambda = 1 / (n_atoms Z sigma_Heitler)
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    explicit Annihilation(const Particle* gamma) : DiscreteProcess("annihil"), gammaDef(gamma) {}

    bool isApplicable(const Particle& particle) const override { return particle.getPdgId() == -11; }

    // [INFO] Public for the validation: Heitler cross section per electron (-> .cpp)
    static double computeCrossSectionPerElectron(double kineticEnergy);

    // In flight
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;

    // At rest: right away (positronium lifetimes, ~0.1-100 ns, are ignored)
    double atRestGetPhysicalInteractionLength(const Track&) override { return 0.0; }
    void atRestDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
