#pragma once

#include "DiscreteProcess.h"

#include <cmath>

// ==============================================================================
// Derived Class: GammaConversion
// Pair production in the field of the nucleus: gamma -> e- e+.
//
// Threshold 2 m_e = 1.022 MeV. Dominant in lead above ~5 MeV; at high energy
// sigma -> (7/9) / (n X0): the mean free path becomes 9/7 X0.
//
// eps = (total energy of the e-) / E_gamma, between m/E and 1 - m/E. Bethe-Heitler
// with Tsai's screening functions (Geant4 Physics Reference Manual):
//   dsigma/deps = alpha r_e^2 Z (Z + xi) { [eps^2 + (1-eps)^2] [Phi1(delta) - F/2]
//                                        + (2/3) eps (1-eps) [Phi2(delta) - F/2] }
//   delta = 136 m_e / (Z^(1/3) E eps (1-eps))   screening: delta -> 0 = complete
//   F = (8/3) ln Z (+ 8 f_c above 50 MeV),  xi(Z): pairs in the field of the atomic electrons
// This formula is used to SAMPLE eps (as G4BetheHeitlerModel does).
//
// The CROSS SECTION comes from the XCOM table (PhotonDataTable), not from the
// integral of this formula. The integral is in the code (computeCrossSectionPerAtom)
// ==============================================================================
class GammaConversion : public DiscreteProcess {
private:
    const Particle* electronDef;
    const Particle* positronDef;

    // Tsai's screening functions
    //   F1 = 3 Phi1 - Phi2,  F2 = 3/2 Phi1 + 1/2 Phi2  (one line each -> here)
    static double screenFunction1(double delta) { return delta > 1.0 ? 42.24 - 8.368 * std::log(delta + 0.952) : 42.392 - delta * (7.796 - 1.961 * delta); }
    static double screenFunction2(double delta) { return delta > 1.0 ? 42.24 - 8.368 * std::log(delta + 0.952) : 41.405 - delta * (5.828 - 0.8945 * delta); }

    // F(Z) of the formula above, with or without the Coulomb correction
    static double computeF(double energy, const Material& material);

protected:
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    GammaConversion(const Particle* electron, const Particle* positron)
        : DiscreteProcess("conv"), electronDef(electron), positronDef(positron) {}

    bool isApplicable(const Particle& particle) const override { return particle.getPdgId() == 22; }

    // The cross section the simulation uses: XCOM, per atom (cm^2)
    static double getCrossSectionPerAtom(double energy, const Material& material);
    // Bethe-Heitler: dsigma/deps, and its integral over eps (NOT used for tracking)
    static double computeDifferentialCrossSection(double energy, double eps, const Material& material);
    static double computeCrossSectionPerAtom(double energy, const Material& material);
    static double sampleEnergyFraction(double energy, const Material& material);

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
