#pragma once

#include "DiscreteProcess.h"

// ==============================================================================
// Derived Class: ComptonScattering
// A photon scatters on an atomic electron, treated as FREE and at rest: Klein-Nishina.
//
// Kinematics (eps = E'/E, k = E/m_e):  1 - cos(theta) = (1 - eps) / (eps k)
// so eps runs from 1/(1 + 2k) (back-scatter) to 1 (forward, nothing happens).
// The electron gets T = E - E'.
//
// Cross section per atom = Z * (Klein-Nishina per electron), exact formula.
// ==============================================================================
class ComptonScattering : public DiscreteProcess {
private:
    const Particle* electronDef; // particle type of the recoil electron

protected:
    // lambda = 1 / (n_atoms Z sigma_KN)
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    explicit ComptonScattering(const Particle* electron) : DiscreteProcess("compt"), electronDef(electron) {}

    bool isApplicable(const Particle& particle) const override { return particle.getPdgId() == 22; }

    // Klein-Nishina total cross section per electron (cm^2)
    static double computeCrossSectionPerElectron(double energy);

    // dsigma/deps per electron, the distribution sampleEnergyFraction draws from
    static double computeDifferentialCrossSection(double energy, double eps);

    // eps = E'/E drawn from Klein-Nishina (G4KleinNishinaCompton::SampleSecondaries)
    static double sampleEnergyFraction(double energy);

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
