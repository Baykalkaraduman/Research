#pragma once

#include "DiscreteProcess.h"


// ==============================================================================
// Derived Class: PhotoElectricEffect
// The photon is absorbed by an atom; an electron of one shell
// leaves with T = E_gamma - B (B: binding energy of that shell).
// ==============================================================================
class PhotoElectricEffect : public DiscreteProcess {
private:
    const Particle* electronDef; // the photoelectron's particle type

    // Direction of the photoelectron (Sauter-Gavrila)
    Vector3D samplePhotoElectronDirection(double electronKineticEnergy, const Vector3D& photonDirection) const;

protected:
    // lambda = 1 / (mu/rho * rho)
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    explicit PhotoElectricEffect(const Particle* electron) : DiscreteProcess("phot"), electronDef(electron) {}

    bool isApplicable(const Particle& particle) const override { return particle.getPdgId() == 22; }

    // Public for the validation: mu/rho in cm^2/g (log-log interpolation)
    double computeMassAttenuation(double energy, const Material& material) const;

    // Binding energy of the shell that absorbs a photon of this energy
    double getBindingEnergy(double energy, const Material& material) const;

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
