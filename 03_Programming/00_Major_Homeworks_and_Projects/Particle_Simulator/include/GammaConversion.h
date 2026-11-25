#pragma once

#include "DiscreteProcess.h"

#include <cmath>
#include <map>
#include <vector>

// ==============================================================================
// Derived Class: GammaConversion
// Pair production in the field of the nucleus: gamma -> e- e+.
//
// Threshold 2 m_e = 1.022 MeV. Dominant in lead above ~5 MeV; at high energy
// sigma -> (7/9) / (n X0): the mean free path becomes 9/7 X0.
//
// eps = (total energy of the e-) / E_gamma, between m_e/E and 1 - m_e/E.
// Bethe-Heitler with Tsai's screening functions (Geant4 Physics Reference Manual):
//   dsigma/deps = alpha r_e^2 Z (Z + xi) { [eps^2 + (1-eps)^2] [Phi1(delta) - F/2]
//                                        + (2/3) eps (1-eps) [Phi2(delta) - F/2] }
//   delta = 136 m_e / (Z^(1/3) E eps (1-eps))   screening: delta -> 0 = complete
//   F = (8/3) ln Z (+ 8 f_c above 50 MeV),  xi(Z): pairs in the field of atomic e-
//
// When LPM is enabled, dsigma/deps is multiplied by Migdal's R_LPM(E, eps) <= 1,
// suppressing symmetric pairs and scaling the high-energy XCOM plateau by
// sigma_LPM(E) / sigma_BH(E) above E_LPM / 50 (~86 GeV in Pb).
// ==============================================================================
class GammaConversion : public DiscreteProcess {
private:
    const Particle* electronDef;
    const Particle* positronDef;
    bool lpmEnabled = false;

    // Cached LPM-suppressed cross section per atom [cm^2] on a ln(E) grid
    struct LpmPairTable {
        double energyMin = 0.0;
        double lnEnergyMin = 0.0;
        double dLnEnergy = 0.0;
        std::vector<double> sigma;
    };
    mutable std::map<const Material*, LpmPairTable> lpmTables;
    static constexpr double tableEnergyMax = 100000.0 * Units::TeV; // 100 PeV
    static constexpr int binsPerDecade = 20;

    [[nodiscard]] const LpmPairTable& getLpmTable(const Material& material) const;
    [[nodiscard]] static double interpolateLpmTable(const LpmPairTable& table, double energy);

    // Tsai's screening functions: F1 = 3 Phi1 - Phi2, F2 = 3/2 Phi1 + 1/2 Phi2
    [[nodiscard]] static double screenFunction1(double delta) {
        return delta > 1.0 ? 42.24 - 8.368 * std::log(delta + 0.952) : 42.392 - delta * (7.796 - 1.961 * delta);
    }
    [[nodiscard]] static double screenFunction2(double delta) {
        return delta > 1.0 ? 42.24 - 8.368 * std::log(delta + 0.952) : 41.405 - delta * (5.828 - 0.8945 * delta);
    }

    // F(Z) of the formula above, with or without the Coulomb correction
    [[nodiscard]] static double computeF(double energy, const Material& material);

protected:
    [[nodiscard]] double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    GammaConversion(const Particle* electron, const Particle* positron)
        : DiscreteProcess("conv"), electronDef(electron), positronDef(positron) {}

    [[nodiscard]] bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 22;}

    void setLPM(bool enable) {lpmEnabled = enable; lpmTables.clear();}
    [[nodiscard]] bool isLPMEnabled() const {return lpmEnabled;}

    // The cross section the simulation uses: XCOM (optionally scaled by LPM), per atom [cm^2]
    [[nodiscard]] static double getCrossSectionPerAtom(double energy, const Material& material, bool enableLPM = false);

    // Bethe-Heitler (optionally with LPM): dsigma/deps and its integral over eps
    [[nodiscard]] static double computeDifferentialCrossSection(double energy, double eps, const Material& material,
                                                                bool enableLPM = false);
    [[nodiscard]] static double computeCrossSectionPerAtom(double energy, const Material& material,
                                                           bool enableLPM = false);
    [[nodiscard]] static double sampleEnergyFraction(double energy, const Material& material, bool enableLPM = false);

    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};