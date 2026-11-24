#pragma once

#include "DiscreteProcess.h"

// ==============================================================================
// Derived Class: GammaConversionToMuons
// [STAGE 7a] New. A photon converts into a MUON pair in the field of a nucleus:
// gamma Z -> mu+ mu- Z. Equivalent to G4GammaConversionToMuons (H. Burkhardt,
// S. Kelner, R. Kokoulin; Geant4 Physics Reference Manual, "Muon pair production
// by photons").
//
// The same diagram as gamma -> e+ e- (stage 3, GammaConversion), with the muon
// mass: the cross section is ~ (m_e / m_mu)^2 = 1/42750 smaller. Threshold
// 2 m_mu = 211 MeV. In a shower it is a rare process, but it is how
// electromagnetic showers make muons (and a background for muon detectors).
//
// Cross section per atom (r_mu = r_e m_e / m_mu, sqrt(e) = 1.6487):
//
//   sigma = (7/9) 4 alpha Z^2 r_mu^2 ln(1 + W_M C_f E_g)
//
//   W_inf = B Z^(-1/3) m_mu / (D_n m_e)       B = 183, D_n = 1.54 A^0.27
//   W_M   = 1 / (4 D_n sqrt(e) m_mu)
//   E_g   = (1 - 4 m_mu / E)^t [ W_sat^s + E^s ]^(1/s),  W_sat = W_inf / W_M,
//           s = -0.88, t = 1.479 + 0.00799 D_n
//   C_f   = 1 + 0.04 ln(1 + E_c / E),          E_c = -18 + 4347 / (B Z^(-1/3))
//           (Geant4's code uses this number in MeV: E_c = 85 MeV for lead, so C_f
//           changes sigma by < 0.4 % above 1 GeV)
//
// At high energy sigma -> (7/9) 4 alpha Z^2 r_mu^2 ln(W_inf): the complete-screening
// limit, the muon analogue of (7/9) / (n X0). Below 5 m_mu the parametrisation is
// replaced by sigma(5 m_mu) y^2, y = (E - 2 m_mu) / (3 m_mu) (as Geant4 does).
//
// Final state (the Geant4 algorithm, -> GammaConversionToMuons.cpp): the energy
// fraction x+ = E(mu+) / E from 1 - 4/3 x+ x- times a logarithm; then the angles
// of both muons from Kelner's approximate distribution in t, psi and rho.
// ==============================================================================
class GammaConversionToMuons : public DiscreteProcess {
private:
    const Particle* muonPlusDef;
    const Particle* muonMinusDef;
    double crossSectionFactor = 1.0;

protected:
    // 1 / (n sigma), with the threshold behaviour
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    GammaConversionToMuons(const Particle* muPlus, const Particle* muMinus)
        : DiscreteProcess("GammaToMuPair"), muonPlusDef(muPlus), muonMinusDef(muMinus) {}

    bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 22;}

    // A rare process: to see enough muon pairs in a test, Geant4 lets one multiply the cross section (the result must then be weighted by 1/f).
    void setCrossSectionFactor(double factor) {crossSectionFactor = factor;}
    double getCrossSectionFactor() const {return crossSectionFactor;}

    // sigma per atom from the formula above, valid for E >= 5 m_mu
    static double computeCrossSectionPerAtom(double photonEnergy, const Material& material);
    // ... and with the threshold factor below 5 m_mu, as used for tracking
    static double getCrossSectionPerAtom(double photonEnergy, const Material& material);

    // x+ = E(mu+) / E_gamma, the Geant4 sampling
    static double sampleEnergyFraction(double photonEnergy, const Material& material);

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
