#pragma once

#include "DiscreteProcess.h"
#include "Units.h"

#include <utility>
#include <vector>
#include <algorithm>
#include <map>

// ==============================================================================
// Enum: FluctuationModel
// How the energy lost in one step is spread around its mean value.
// ==============================================================================
enum class FluctuationModel {
    none,       // continuous energy losing
    gauss,      // Thick layers, where many collisions add up
    landau      // Thin layers, where a few hard collisions happens
};

// ==============================================================================
// Derived Class: Ionization
// Mocks the continuous energy loss of a charged particle in a medium.
// ==============================================================================
class Ionization : public DiscreteProcess{
private:
    // ==========================================================================
    // Step limitation
    // ==========================================================================
    // dE/dx changes with energy, so one step must not lose too much energy, or
    // using the pre-step dE/dx for the whole step becomes inaccurate.
    double dRmaxRange = 0.2;                        // max fraction of the range per step
    double finalRange = 0.1 * Units::cm;            // below this range: allowed to stop in one step
    double lowestKineticEnergy = 1.0 * Units::keV;  // tracking cut: below this, deposit everything and stop
    double linLossLimit = 0.01;                     // limit for dE/dx * step
    FluctuationModel fluctuationModel = FluctuationModel::gauss; // Fluctuation Model, predefined with gaussian

    // knock-on electrons that have enough kinetic energy to cause ionization away from the primary particle
    const Particle* deltaElectronDef = nullptr;
    double deltaRayCut = GeoDefs::kInfinity;

    // the largest sigma over the coming step, and how often sigma at the end of a step came out larger than it (should be ~0)
    double sigmaMax = 0.0;
    long sigmaMaxViolations = 0;

    // Energy of the (only) maximum of sigma(T) per (particle, material), found once by a scan; kInfinity if sigma only rises
    std::map<std::pair<const Particle*, const Material*>, double> sigmaMaxEnergies;
    double getSigmaMaxEnergies(const Particle& particle, const Material& material);

    // ==========================================================================
    // Range table
    // ==========================================================================
    // R(T) on a logarithmic energy grid. ln(T) is equally spaced, so the bin of any T is found with one division (no search).
    struct RangeTable {
        double lnEnergyMin = 0.0;         // ln of the first grid energy
        double dLnEnergy = 0.0;           // grid step in ln(T)
        std::vector<double> range;        // R at each grid energy (cm)
        std::vector<double> slope;        // dR/d(lnT) = T/S at each grid energy (cm)
    };
    // One table per (particle type, material) pair. Pointers are fine as keys:
    // Particle and Material objects live for the whole run.
    std::map<std::pair<const Particle*, const Material*>, RangeTable> rangeTables;

    static constexpr double tableEnergyMin = 1.0 * Units::keV;
    static constexpr double tableEnergyMax = 100.0 * Units::TeV;
    static constexpr int binsPerDecade = 50;

    // Builds the table on first use, then returns it
    const RangeTable& getRangeTable(const Particle& particle, const Material& material);

    // ==========================================================================
    // Stopping-power models
    // ==========================================================================
    // Heavy charged particles (mu, tau, pi, K, p): Bethe-Bloch, PDG eq. 34.5
    double computeBetheBlochDEDX(const Particle& particle, double kineticEnergy, const Material& material) const;

    // Electrons and positrons: Berger-Seltzer (ICRU 37), collision part only
    double computeElectronDEDX(const Particle& particle, double kineticEnergy, const Material& material) const;

    // Below this energy the models are not valid; dE/dx is then scaled with sqrt(T)
    double getLowEnergyLimit(const Particle& particle) const;

    // e- or e+ ? (They need a different formula and a different Tmax.)
    static bool isElectronOrPositron(const Particle& particle) { return particle.getPdgId() == 11 || particle.getPdgId() == -11; }

    // The stopping power of the CONTINUOUS part: restricted if delta rays are on, the full one otherwise.
    // The range table, the step limit and the along-step loss all use this one.
    double computeContinuousDEDX(const Particle& particle, double kineticEnergy, const Material& material) const {return deltaElectronDef != nullptr ? computeRestrictedDEDX(particle, kineticEnergy, material, deltaRayCut) : computeDEDX(particle, kineticEnergy, material);}

    // Delta-ray helpers
    // Largest delta-ray energy: Tmax, but never more than T
    double getMaxDeltaRayEnergy(const Particle& particle, double kineticEnergy) const { return std::min(computeMaxSecondaryEnergy(particle, kineticEnergy), kineticEnergy); }
    // P = (K/2) z^2 (Z/A) rho / beta^2, energy / length
    static double computeDeltaRayPrefactor(const Particle& particle, double kineticEnergy, const Material& material);
    // The shape z(T') of dsigma/dT' = P z / T'^2 (dimensionless, ~1)
    double computeDeltaRayShape(const Particle& particle, double kineticEnergy, double deltaEnergy) const;


    // ==========================================================================
    // Landau sampling helpers
    // ==========================================================================
    double sampleLandauLoss(double meanLoss, double xi, double lambdaMax) const;

protected:
    // Mean free path used to sample the next delta ray: 1 / sigmaMax, the largest sigma(T_cut) over the coming step (integral approach)
    double getMeanFreePath(const Track&, const Material&) const override { return sigmaMax > 0.0 ? 1.0 / sigmaMax : GeoDefs::kInfinity; }

public:
    // Constructor
    Ionization() : DiscreteProcess("Ionization") {}

    // Delta rays above cut, as the electron tracks. The range table depends on the cut, so it is rebuilt. Set before the run (like Bremsstrahlung).
    void setDeltaRayProduction(const Particle* electron, double cut) {deltaElectronDef = electron; deltaRayCut = cut; rangeTables.clear(); sigmaMaxEnergies.clear();}
    double getDeltaRayCut() const {return deltaRayCut;}
    long getSigmaMaxViolations() const {return sigmaMaxViolations;}

    bool isApplicable(const Particle& particle) const override { return particle.getCharge() != 0.0; }

    // Choose the fluctuation model
    void setFluctuationModel(FluctuationModel model) {fluctuationModel = model;}
    FluctuationModel getFluctuationModel() const {return fluctuationModel;}

    // ==========================================================================
    // Physics quantities. Public, so the validation program can compare them with the PDG tables directly
    // ==========================================================================
    // Mean (unrestricted) stopping power in GeV/cm
    double computeDEDX(const Particle& particle, double kineticEnergy, const Material& material) const;

    // Largest energy one collision can give to an electron
    double computeMaxSecondaryEnergy(const Particle& particle, double kineticEnergy) const;

    // Restricted stopping power: only collisions with T' < cut
    double computeRestrictedDEDX(const Particle& particle, double kineticEnergy, const Material& material, double cut) const;

    double computeDeltaRayDifferentialCrossSection(const Particle& particle, double kineticEnergy, double deltaEnergy,
                                                   const Material& material) const;
    // sigma(T' > cut) in 1 / length
    double computeDeltaRayCrossSection(const Particle& particle, double kineticEnergy, const Material& material, double cut) const;
    // integral_cut^Tmax T' dsigma/dT', energy / length: what delta rays carries
    double computeDeltaRayEnergyLoss(const Particle& particle, double kineticEnergy, const Material& material, double cut) const;
    // T' drawn from dsigma/dT' between cut and Tmax
    double sampleDeltaRayEnergy(const Particle& particle, double kineticEnergy, double cut) const;

    // CSDA range in cm from the table
    double getRange(const Particle& particle, double kineticEnergy, const Material& material);

    // Inverse of getRange: the kinetic energy that has exactly this range
    double getKineticEnergyFromRange(const Particle& particle, double range, const Material& material);

    // Actual loss of one step, drawn around meanLoss
    double sampleFluctuations(const Particle& particle, double kineticEnergy, const Material& material,
                              double meanLoss, double stepLength) const;

    // Landau most probable loss Delta_p for a layer of thickness stepLength (PDG eq. 34.12)
    double computeMostProbableLoss(const Particle& particle, double kineticEnergy, const Material& material,
                                   double stepLength) const;

    double alongStepGetPhysicalInteractionLength(const Track& track, const Material& material) override;
    void alongStepDoIt(Track& track, const Material& material, double stepLength) override;

    // Post-step half: distance to the next delta ray. Without delta rays:
    // kInfinity at once (and NO random number is drawn, so programs that do not switch delta rays on keep their random sequence).
    double postStepGetPhysicalInteractionLength(const Track& track, const Material& material) override;

    // One delta ray
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;

};
