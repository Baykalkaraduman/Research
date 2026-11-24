#pragma once

#include "DiscreteProcess.h"

class AtomicDeexcitation;

// ==============================================================================
// Derived Class: PhotoElectricEffect
// Absorbs a photon on a sampled atomic subshell (ENDF/EPDL), emits a Sauter-
// Gavrila photoelectron with T = E_gamma - B_shell, and triggers atomic
// relaxation (fluorescence / Auger cascade) when an AtomicDeexcitation is set.
// ==============================================================================
class PhotoElectricEffect : public DiscreteProcess {
private:
    // --------------------------------------------------------------------------
    // Member Data & Internal Helpers
    // --------------------------------------------------------------------------
    const Particle* electronDef;
    const AtomicDeexcitation* deexcitation = nullptr;

    [[nodiscard]] Vector3D samplePhotoElectronDirection(double electronKineticEnergy,
                                                          const Vector3D& photonDirection) const;

protected:
    // --------------------------------------------------------------------------
    // DiscreteProcess Hooks
    // --------------------------------------------------------------------------
    [[nodiscard]] double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    // --------------------------------------------------------------------------
    // Construction & Configuration
    // --------------------------------------------------------------------------
    explicit PhotoElectricEffect(const Particle* electron) : DiscreteProcess("phot"), electronDef(electron) {}

    void setAtomicDeexcitation(const AtomicDeexcitation* deex) {deexcitation = deex;}

    // --------------------------------------------------------------------------
    // Process & Physics Interface
    // --------------------------------------------------------------------------
    [[nodiscard]] bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 22;}

    [[nodiscard]] double computeMassAttenuation(double energy, const Material& material) const;
    [[nodiscard]] double getBindingEnergy(double energy, const Material& material) const;

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};