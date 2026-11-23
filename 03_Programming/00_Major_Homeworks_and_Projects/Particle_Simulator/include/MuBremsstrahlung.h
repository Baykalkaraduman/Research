#pragma once

#include "EnergyLossProcess.h"

#include <cstdlib>

// ==============================================================================
// Derived Class: MuBremsstrahlung
// [STAGE 7a] A muon is deflected in the field of a nucleus (or atomic electron)
// and radiates a photon of energy k = eps.
// Equivalent to G4MuBremsstrahlung with G4MuBremsstrahlungModel.
//
// Same physics as electron bremsstrahlung (stage 3), but the radiating particle
// is 207 times heavier: the cross section scales as (m_e / M)^2 = 1/42750.
// Consequently, a muon in lead loses energy mostly by ionization up to ~140 GeV
// (muon critical energy), compared to ~7 MeV for an electron.
//
// Differential cross section (Kelner, Kokoulin, Petrukhin 1995; Geant4 Physics
// Reference Manual, "Muon bremsstrahlung"), v = k / E, E = total energy:
//
//   dsigma/dk = (16/3) alpha r_mu^2 Z (1 - v + 3/4 v^2) (Z Phi_n + Phi_e) / k,
//   r_mu = r_e m_e / M   (classical radius of the muon)
//
// With minimum momentum transfer delta = M^2 v / (2 (E - k)) and sqrt(e) = 1.6487:
//
//   Phi_n = ln[ B Z^(-1/3) (M + delta (D_n' sqrt(e) - 2)) / (D_n' (m_e + delta sqrt(e) B Z^(-1/3))) ]
//           Nucleus: atomic screening (B = 183) and finite nuclear size
//           (D_n = 1.54 A^0.27, D_n' = D_n^(1 - 1/Z)).
//
//   Phi_e = ln[ B' Z^(-2/3) M / ((1 + delta M / (m_e^2 sqrt(e))) (m_e + delta sqrt(e) B' Z^(-2/3))) ]
//           Atomic electrons (B' = 1429), valid only below k_max' = E / (1 + M^2 / (2 m_e E)).
//
// Both logarithms are clamped at 0 from below. Without screening and nuclear
// size (large delta), Phi_n -> ln(M / delta): the classic Bethe-Heitler result.
//
// [INFO] Production cut: photons below k_c = 100 keV are treated as continuous
//        loss (Geant4's default 0.7 mm range cut gives ~100 keV in lead). Since
//        k * dsigma/dk is nearly flat, this is only ~k_c / E of the brems loss
//        (~1e-6 at 100 GeV).
// ==============================================================================
class MuBremsstrahlung : public EnergyLossProcess {
private:
    // --------------------------------------------------------------------------
    // Constants & Member Data
    // --------------------------------------------------------------------------
    // G4MuBremsstrahlungModel defaults: minThreshold (0.9 keV) and lowestKinEnergy (0.1 GeV)
    static constexpr double minThreshold = 0.9 * Units::keV;
    static constexpr double lowestKineticEnergy = 0.1 * Units::GeV;

    const Particle* gammaDef;

protected:
    // --------------------------------------------------------------------------
    // EnergyLossProcess Hooks
    // --------------------------------------------------------------------------
    double getMinEnergyTransfer(const Particle&, const Material&) const override {return minThreshold;}
    double getMaxEnergyTransfer(const Particle&, double kineticEnergy, const Material&) const override {return kineticEnergy;}
    double getLowestKineticEnergy(const Particle&) const override {return lowestKineticEnergy;}

    // Splits numerical integration at k_max' where Phi_e drops to 0
    std::vector<double> getIntegrationBreakpoints(const Particle& particle, double kineticEnergy, const Material& material) const override;

    // Produces one radiated photon (-> MuBremsstrahlung.cpp)
    void produceSecondaries(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries, double transfer) override;

public:
    // --------------------------------------------------------------------------
    // Construction
    // --------------------------------------------------------------------------
    explicit MuBremsstrahlung(const Particle* gamma) : EnergyLossProcess("muBrems", 100.0 * Units::keV), gammaDef(gamma) {}

    // --------------------------------------------------------------------------
    // Process & Physics Interface
    // --------------------------------------------------------------------------
    // Applicable to mu+ and mu- (PDG ID +/-13)
    bool isApplicable(const Particle& particle) const override {return std::abs(particle.getPdgId()) == 13;}

    // Kelner-Kokoulin-Petrukhin differential cross section per atom (-> MuBremsstrahlung.cpp)
    double computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double transfer, const Material& material) const override;
};