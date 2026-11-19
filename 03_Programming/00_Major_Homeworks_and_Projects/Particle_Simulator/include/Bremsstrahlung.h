#pragma once

#include "DiscreteProcess.h"

#include <map>
#include <vector>

// ==============================================================================
// Derived Class: Bremsstrahlung
// An e- or e+ is deflected by the nuclear (and atomic-electron)
// field and radiates a photon of energy k.
//
// In lead this is THE energy loss of an electron above ~10 MeV (critical energy),
// and with pair production it builds the electromagnetic shower.
//
// dsigma/dk ~ 1/k: infinitely many soft photons. So, like Geant4, the photons are
// split at a production cut k_c:
//   k <  k_c : not produced; their energy is a CONTINUOUS loss (along-step), left
//              along the track (it would be absorbed close by anyway)
//   k >= k_c : produced as real photons (post-step, discrete), cross section
//              sigma(T, k_c) = integral_{k_c}^{T} dsigma/dk dk
//
// Differential cross section (Tsai, complete + partial screening; y = k/E, E total):
//   dsigma/dk = (16/3) alpha r_e^2 Z^2 / k * { (1 - y + 3/4 y^2) [ Phi1/4 - ln Z/3 - f_c
//                                                              + (Psi1/4 - 2/3 ln Z) / Z ]
//                                             + (1 - y)/8 [ Phi1m2 + Psi1m2 / Z ] }
// ==============================================================================
class Bremsstrahlung : public DiscreteProcess {
private:
    const Particle* gammaDef;
    double productionCut = 100.0 * Units::keV; // k_c

    // Two tables per material on a ln(T) grid: sigma(k > k_c) per atom and the
    // continuous loss of k < k_c, per atom (energy * cm^2)
    struct BremsTable {
        double lnEnergyMin = 0.0;
        double dLnEnergy = 0.0;
        std::vector<double> sigma;
        std::vector<double> softLoss;
    };
    mutable std::map<const Material*, BremsTable> tables;
    static constexpr double tableEnergyMin = 1.0 * Units::keV;
    static constexpr double tableEnergyMax = 100.0 * Units::TeV;
    static constexpr int binsPerDecade = 20;
    const BremsTable& getTable(const Material& material) const;
    static double interpolate(const BremsTable& table, const std::vector<double>& values, double energy);

    // The curly bracket above, dimensionless
    static double computeDxSection(double kineticEnergy, double photonEnergy, const Material& material);

protected:
    double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    explicit Bremsstrahlung(const Particle* gamma) : DiscreteProcess("eBrem"), gammaDef(gamma) {}

    bool isApplicable(const Particle& particle) const override { return particle.getPdgId() == 11 || particle.getPdgId() == -11; }

    // The cut must be set before the run: the tables depend on it
    void setProductionCut(double cut) { productionCut = cut; tables.clear(); }
    double getProductionCut() const { return productionCut; }

    // dsigma/dk per atom (cm^2 / energy)
    static double computeDifferentialCrossSection(double kineticEnergy, double photonEnergy, const Material& material);
    // sigma per atom for photons with k > kMin
    static double computeCrossSectionPerAtom(double kineticEnergy, double kMin, const Material& material);
    // integral_0^kMax k dsigma/dk, per atom; with kMax = T it gives the radiative
    // stopping power (ESTAR's "radiative" column) times 1/n
    static double computeEnergyLossPerAtom(double kineticEnergy, double kMax, const Material& material);
    // Full radiative stopping power, energy / cm
    static double computeRadiativeDEDX(double kineticEnergy, const Material& material);
    // k drawn from dsigma/dk between k_c and T
    static double samplePhotonEnergy(double kineticEnergy, double kMin, const Material& material);

    // Continuous part: loss of the photons below k_c
    void alongStepDoIt(Track& track, const Material& material, double stepLength) override;

    // Discrete part: emit one photon above k_c
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};
