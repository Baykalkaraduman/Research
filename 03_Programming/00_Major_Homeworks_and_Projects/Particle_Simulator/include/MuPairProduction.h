#pragma once

#include "EnergyLossProcess.h"

#include <algorithm>
#include <cstdlib>

// ==============================================================================
// Derived Class: MuPairProduction
// [STAGE 7a] New. Direct e+e- pair production by a muon in the field of a
// nucleus: mu Z -> mu e+ e- Z. The muon gives eps = E(e+) + E(e-) (total
// energies) to the pair.
// Equivalent to G4MuPairProduction with G4MuPairProductionModel.
//
// Think of it as the muon's virtual photon converting into a pair (the
// "trident" process). It is the LARGEST energy loss of a muon in lead above
// ~100 GeV, but it is made of many small transfers: eps dsigma/deps peaks at
// eps ~ 1e-3 E. Brems is fewer, harder transfers.
//
// Cross section (Kokoulin & Petrukhin 1971; R.P. Kokoulin's formula as coded
// in G4MuPairProductionModel), v = eps / E:
//
//   dsigma/dv = (4/(3 pi)) (alpha r_e)^2 Z (Z + zeta) (1 - v)/v
//               * integral_0^{rho_max} [ Phi_e + (m_e/M)^2 Phi_mu ] drho
//
//   rho = (E+ - E-) / eps   the energy asymmetry of the pair
//   rho_max = [1 - 6 M^2 / (E^2 (1 - v))] sqrt(1 - 4 m_e / eps)
//   Phi_e, Phi_mu: the diagrams where the e- or the muon line takes the
//                  momentum transfer from the nucleus, with screening (B = 183)
//                  and nuclear size; long formulas, see MuPairProduction.cpp
//   zeta: pairs in the field of the atomic electrons (like Z(Z+1) elsewhere)
//
// The rho integral is done with 8-point Gauss in ln(1 - rho) (as Geant4 does),
// so ONE value of dsigma/deps costs 8 evaluations; that is why the tables matter.
//
// Limits: 4 m_e <= eps <= T + M (1 - (3/4) sqrt(e) Z^(1/3)). Lowest muon energy
// max(0.85 GeV, 8 M) = 0.85 GeV
// ==============================================================================
class MuPairProduction : public EnergyLossProcess {
private:
    const Particle* electronDef;
    const Particle* positronDef;

    static constexpr double lowestKineticEnergy = 0.85 * Units::GeV;

protected:
    // 4 m_e
    double getMinEnergyTransfer(const Particle& particle, const Material& material) const override;
    // T + M (1 - 0.75 sqrt(e) Z^(1/3))
    double getMaxEnergyTransfer(const Particle& particle, double kineticEnergy, const Material& material) const override;
    double getLowestKineticEnergy(const Particle& particle) const override {return std::max(lowestKineticEnergy, 8.0 * particle.getMass());}

    // e- and e+
    void produceSecondaries(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries,
                            double transfer) override;

public:
    MuPairProduction(const Particle* electron, const Particle* positron)
    : EnergyLossProcess("muPairProd", 1.0 * Units::MeV), electronDef(electron), positronDef(positron) {}

    bool isApplicable(const Particle& particle) const override {return std::abs(particle.getPdgId()) == 13;}

    // Kokoulin's formula
    double computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double transfer, const Material& material) const override;
};
