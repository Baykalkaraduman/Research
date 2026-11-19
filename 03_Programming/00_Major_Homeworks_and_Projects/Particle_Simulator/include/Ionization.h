#pragma once

#include "PhysicsProcess.h"
#include "Units.h"

#include <cmath>
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
class Ionization : public PhysicsProcess {
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

    // ==========================================================================
    // Range table
    // ==========================================================================
    // R(T) on a logarithmic energy grid. ln(T) is equally spaced, so the bin of any T is found with one division (no search).
    struct RangeTable {
        double lnEnergyMin = 0.0;         // ln of the first grid energy
        double dLnEnergy = 0.0;           // grid step in ln(T)
        std::vector<double> range;        // R at each grid energy (cm)
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

    // ==========================================================================
    // Landau sampling helpers
    // ==========================================================================
    double sampleLandauLoss(double meanLoss, double xi, double lambdaMax) const;

public:
    // Constructor
    Ionization() : PhysicsProcess("Ionization") {}

    bool isApplicable(const Particle& particle) const override { return particle.getCharge() != 0.0; }

    // Choose the fluctuation model
    void setFluctuationModel(FluctuationModel model) {fluctuationModel = model;}
    FluctuationModel getFluctuationModel() const {return fluctuationModel;}

    // ==========================================================================
    // Physics quantities. Public, so the validation program can compare them with the PDG tables directly
    // ==========================================================================
    // Mean (unrestricted) stopping power in GeV/cm (G4VEmModel::ComputeDEDX)
    double computeDEDX(const Particle& particle, double kineticEnergy, const Material& material) const;

    // Largest energy one collision can give to an electron
    double computeMaxSecondaryEnergy(const Particle& particle, double kineticEnergy) const;

    // CSDA range in cm from the table (G4VEnergyLossProcess::GetRange)
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
};
