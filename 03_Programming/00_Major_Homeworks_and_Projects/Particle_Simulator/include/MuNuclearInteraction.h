#pragma once

#include "EnergyLossProcess.h"
#include "PhysicalConstants.h"

#include <cstdlib>

// ==============================================================================
// Derived Class: MuNuclearInteraction
// [STAGE 7c] A muon interacts inelastically with a nucleus via a virtual photon
// (photonuclear / deep-inelastic scattering), transferring energy eps = nu and
// producing secondary hadrons (nucleon + pion(s)).
//
// Differential cross section (Bezrukov & Bugaev 1981; Geant4 Physics Reference
// Manual, "Muon Photonuclear Interaction"), v = eps / E, E = total energy:
//
//   dsigma/deps = (alpha / 2pi) A sigma_gammaN(eps) (v^2 / eps)
//                 * { 0.75 G(x) [ kappa ln(1 + m1^2 / t) - kappa m1^2 / (m1^2 + t) - 2 M^2 / t ]
//                   + 0.25 [ kappa ln(1 + m2^2 / t) - 2 M^2 / t ]
//                   + (M^2 / 2t) [ 0.75 G(x) m1^2 / (m1^2 + t) + 0.25 (m2^2 / t) ln(1 + t / m2^2) ] }
//
// with:
//   kappa = 1 - 2/v + 2/v^2,   t = M^2 v^2 / (1 - v)
//   sigma_gammaN(eps) = [ 114.3 + 1.647 ln^2(0.0213 eps_GeV) ] ub
//   G(x) = (3 / x^3) (x^2 / 2 - 1 + (1 + x) e^-x),   x = 0.00282 A^(1/3) sigma_gammaN[ub]
//
// Transfer limits:
//   eps_min = 0.2 GeV  (above single-pion photoproduction threshold)
//   eps_max = min(T, E - (m_p / 2) (1 + M^2 / m_p^2))
//
// [INFO] Production cut: transfers below eps_c = 1 GeV are treated as continuous
//        along-step energy loss; transfers above eps_c produce explicit hadrons.
// ==============================================================================
class MuNuclearInteraction : public EnergyLossProcess {
private:
    // --------------------------------------------------------------------------
    // Constants & Member Data
    // --------------------------------------------------------------------------
    static constexpr double minThreshold = 0.2 * Units::GeV;
    static constexpr double lowestKineticEnergy = 1.0 * Units::GeV;

    const Particle* protonDef;
    const Particle* piPlusDef;
    const Particle* piMinusDef;
    const Particle* pi0Def;

protected:
    // --------------------------------------------------------------------------
    // EnergyLossProcess Hooks
    // --------------------------------------------------------------------------
    double getMinEnergyTransfer(const Particle&, const Material&) const override {return minThreshold;}
    double getMaxEnergyTransfer(const Particle& particle, double kineticEnergy, const Material& material) const override;
    double getLowestKineticEnergy(const Particle&) const override {return lowestKineticEnergy;}

    // Samples virtual-photon Q^2, updates scattered muon, and boosts secondary hadrons
    void produceSecondaries(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries, double transfer) override;

public:
    // --------------------------------------------------------------------------
    // Construction
    // --------------------------------------------------------------------------
    MuNuclearInteraction(const Particle* proton, const Particle* piPlus, const Particle* piMinus, const Particle* pi0 = nullptr)
        : EnergyLossProcess("muNucl", 1.0 * Units::GeV), protonDef(proton), piPlusDef(piPlus), piMinusDef(piMinus), pi0Def(pi0) {}

    // --------------------------------------------------------------------------
    // Process & Physics Interface
    // --------------------------------------------------------------------------
    // Applicable to mu+ and mu- (PDG ID +/-13)
    bool isApplicable(const Particle& particle) const override {return std::abs(particle.getPdgId()) == 13;}

    // Bezrukov-Bugaev differential cross section per atom (-> MuNuclearInteraction.cpp)
    double computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double transfer, const Material& material) const override;
};