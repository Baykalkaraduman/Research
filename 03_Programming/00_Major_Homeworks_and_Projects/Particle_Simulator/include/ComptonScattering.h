#pragma once

#include "DiscreteProcess.h"

// ==============================================================================
// Derived Class: ComptonScattering
// Incoherent photon scattering on atomic electrons.
// Cross section: ENDF/B-VIII.0 MF=23 MT=504 (via PhotoElectricShellData).
// Kinematics: Klein-Nishina rejection modulated by S(x)/Z (MF=27 MT=504),
// with recoil electron kinetic energy T = E - E'.
// Free-electron static methods are kept for validateEM (Klein-Nishina checks).
// ==============================================================================
class ComptonScattering : public DiscreteProcess {
private:
    // --------------------------------------------------------------------------
    // Member Data
    // --------------------------------------------------------------------------
    const Particle* electronDef;

protected:
    // --------------------------------------------------------------------------
    // DiscreteProcess Hooks
    // --------------------------------------------------------------------------
    [[nodiscard]] double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    // --------------------------------------------------------------------------
    // Construction & Process Interface
    // --------------------------------------------------------------------------
    explicit ComptonScattering(const Particle* electron) : DiscreteProcess("compt"), electronDef(electron) {}

    [[nodiscard]] bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 22;}

    // --------------------------------------------------------------------------
    // Cross Sections & Energy-Fraction Sampling (public for validation)
    // --------------------------------------------------------------------------
    [[nodiscard]] static double computeCrossSectionPerElectron(double energy);
    [[nodiscard]] static double computeDifferentialCrossSection(double energy, double eps);
    [[nodiscard]] double computeCrossSectionPerAtom(double energy, const Material& material) const;

    [[nodiscard]] static double sampleEnergyFraction(double energy);
    [[nodiscard]] static double sampleEnergyFraction(double energy, int z);

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};