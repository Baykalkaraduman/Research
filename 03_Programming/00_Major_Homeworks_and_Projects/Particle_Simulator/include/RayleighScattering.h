#pragma once

#include "DiscreteProcess.h"

// ==============================================================================
// Derived Class: RayleighScattering
// Coherent (elastic) photon scattering on bound atomic electrons.
// Cross section: ENDF/B-VIII.0 MF=23 MT=502.
// Angular distribution: dsigma/dOmega ~ (1 + cos^2(theta))/2 * F^2(x),
// sampled via cumulative F^2(x) in x^2 (MF=27 MT=502) with rejection on
// (1 + cos^2(theta))/2. Photon energy is unchanged and no secondaries are emitted.
// ==============================================================================
class RayleighScattering : public DiscreteProcess {
protected:
    // --------------------------------------------------------------------------
    // DiscreteProcess Hooks
    // --------------------------------------------------------------------------
    [[nodiscard]] double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    // --------------------------------------------------------------------------
    // Construction & Process Interface
    // --------------------------------------------------------------------------
    RayleighScattering() : DiscreteProcess("Rayl") {}

    [[nodiscard]] bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 22;}

    // --------------------------------------------------------------------------
    // Physics Queries & Sampling (public for validation)
    // --------------------------------------------------------------------------
    [[nodiscard]] double computeCrossSectionPerAtom(double photonEnergy, const Material& material) const;
    [[nodiscard]] double computeMassAttenuation(double photonEnergy, const Material& material) const;
    [[nodiscard]] double sampleCosTheta(double photonEnergy, int z) const;

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};