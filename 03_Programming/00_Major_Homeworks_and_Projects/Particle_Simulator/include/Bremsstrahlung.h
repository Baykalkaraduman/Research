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
//
// High-energy medium suppression:
//   - setLPM(bool) and setDielectricSuppression(bool) operate independently.
//   - Both default to false so ValidationEMApp reproduces NIST ESTAR directly.
// ==============================================================================
class Bremsstrahlung : public DiscreteProcess {
private:
    const Particle* gammaDef;
    double productionCut = 100.0 * Units::keV; // k_c
    bool lpmEnabled = false;
    bool dielectricEnabled = false;

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
    static constexpr double tableEnergyMax = 10000.0 * Units::TeV; // 10 PeV for LPM regime
    static constexpr int binsPerDecade = 20;

    [[nodiscard]] const BremsTable& getTable(const Material& material) const;
    [[nodiscard]] static double interpolate(const BremsTable& table, const std::vector<double>& values, double energy);

    // The curly bracket above (dimensionless), optionally modulated by LPM and/or dielectric suppression
    [[nodiscard]] static double computeDxSection(double kineticEnergy, double photonEnergy, const Material& material,
                                                 bool enableLPM = false, bool includeDielectric = false);

protected:
    [[nodiscard]] double getMeanFreePath(const Track& track, const Material& material) const override;

public:
    explicit Bremsstrahlung(const Particle* gamma) : DiscreteProcess("eBrem"), gammaDef(gamma) {}

    [[nodiscard]] bool isApplicable(const Particle& particle) const override {
        return particle.getPdgId() == 11 || particle.getPdgId() == -11;
    }

    // Configuration (tables depend on cut and suppression flags, so setters invalidate cache)
    void setProductionCut(double cut) {productionCut = cut; tables.clear();}
    [[nodiscard]] double getProductionCut() const {return productionCut;}

    void setLPM(bool enable) {lpmEnabled = enable; tables.clear();}
    [[nodiscard]] bool isLPMEnabled() const {return lpmEnabled;}

    void setDielectricSuppression(bool enable) {dielectricEnabled = enable; tables.clear();}
    [[nodiscard]] bool isDielectricSuppressionEnabled() const {return dielectricEnabled;}

    // dsigma/dk per atom (cm^2 / energy)
    [[nodiscard]] static double computeDifferentialCrossSection(double kineticEnergy, double photonEnergy,
                                                                const Material& material, bool enableLPM = false,
                                                                bool includeDielectric = false);
    // sigma per atom for photons with k > kMin
    [[nodiscard]] static double computeCrossSectionPerAtom(double kineticEnergy, double kMin, const Material& material,
                                                           bool enableLPM = false, bool includeDielectric = false);
    // integral_0^kMax k dsigma/dk, per atom; with kMax = T it gives the radiative
    // stopping power (ESTAR's "radiative" column) times 1/n
    [[nodiscard]] static double computeEnergyLossPerAtom(double kineticEnergy, double kMax, const Material& material,
                                                         bool enableLPM = false, bool includeDielectric = false);
    // Full radiative stopping power, energy / cm
    [[nodiscard]] static double computeRadiativeDEDX(double kineticEnergy, const Material& material,
                                                     bool enableLPM = false, bool includeDielectric = false);
    // k drawn from dsigma/dk between kMin and T
    [[nodiscard]] static double samplePhotonEnergy(double kineticEnergy, double kMin, const Material& material,
                                                   bool enableLPM = false, bool includeDielectric = false);

    // Continuous part: loss of the photons below k_c
    void alongStepDoIt(Track& track, const Material& material, double stepLength) override;

    // Discrete part: emit one photon above k_c
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;
};