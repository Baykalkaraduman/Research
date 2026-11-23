#pragma once

#include "DiscreteProcess.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

// ==============================================================================
// Abstract Class: EnergyLossProcess
// [STAGE 7a] Common machinery for the three radiative energy losses of a muon:
// bremsstrahlung, e+e- pair production, and photonuclear interaction.
// Equivalent to G4VEnergyLossProcess (non-ionization part).
//
// Transfers (eps) are split at a production cut (eps_c):
//   eps <  eps_c : CONTINUOUS (along-step)
//                  Restricted loss dE/dx_restr = n * int_{eps_min}^{eps_c} eps * (dsigma/deps) deps
//                  is subtracted over the step; no secondaries are produced.
//   eps >= eps_c : DISCRETE (post-step)
//                  With sigma(eps_c) = int_{eps_c}^{eps_max} (dsigma/deps) deps per atom,
//                  eps is sampled from dsigma/deps and secondaries are created by the derived class.
//
// Restricted + discrete always equals the full loss n * int eps dsigma;
// the cut only controls where the energy goes.
//
// Derived classes define dsigma/deps, the transfer limits, and secondary
// production. Integration, lookup tables, sampling, and stepping live here.
// ==============================================================================
class EnergyLossProcess : public DiscreteProcess {
private:
    // --------------------------------------------------------------------------
    // Constants & Internal Types
    // --------------------------------------------------------------------------
    static constexpr double tableEnergyMax = 1.0e6 * Units::GeV; // 1 PeV, matching PDG tables
    static constexpr int binsPerDecade = 20;

    // Tabulated on a ln(T) grid from getLowestKineticEnergy() to tableEnergyMax.
    struct LossTable {
        double lnEnergyMin = 0.0;
        double dLnEnergy = 0.0;
        std::vector<double> sigma;          // n * sigma(eps > eps_c), 1/cm
        std::vector<double> restrictedDEDX; // n * int_{eps < eps_c} eps dsigma, energy/cm
        std::vector<double> majorant;       // Upper bound of eps * dsigma/deps on [eps_c, eps_max] per atom
    };

    using TableKey = std::pair<const Particle*, const Material*>;

    // --------------------------------------------------------------------------
    // Internal Helpers
    // --------------------------------------------------------------------------
    // Builds and caches the table on first use.
    const LossTable& getTable(const Particle& particle, const Material& material) const;

    // Linear interpolation in ln(T); clamped to the last node above tableEnergyMax.
    static double interpolate(const LossTable& table, const std::vector<double>& values, double energy);

    // Largest eps * dsigma/deps on [lo, hi], scaled by a safety factor.
    double computeMajorant(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const;

    // Samples eps from dsigma/deps on [lo, hi] using rejection sampling with majorant fMax.
    double sampleWithMajorant(const Particle& particle, double kineticEnergy, double lo, double hi, double fMax, const Material& material) const;

    // --------------------------------------------------------------------------
    // Member Data
    // --------------------------------------------------------------------------
    // Configuration
    double productionCut = 0.0; // eps_c

    // Stepping state (saved by alongStepGPIL for the trapezoid rule in alongStepDoIt)
    double preStepDEDX = 0.0;

    // Cached physics tables per (particle, material)
    mutable std::map<TableKey, LossTable> tables;

    // Diagnostic counters (both should remain 0 in valid runs)
    mutable long majorantViolations = 0;
    long sigmaMaxViolations = 0;

protected:
    // --------------------------------------------------------------------------
    // Derived-Class Physics Hooks
    // --------------------------------------------------------------------------
    // Smallest possible transfer (e.g., 4 m_e for pair production) and largest at kineticEnergy.
    virtual double getMinEnergyTransfer(const Particle& particle, const Material& material) const = 0;

    virtual double getMaxEnergyTransfer(const Particle& particle, double kineticEnergy, const Material& material) const = 0;

    // Below this kinetic energy, the process is switched off (model invalid).
    virtual double getLowestKineticEnergy(const Particle& particle) const = 0;

    // Executes one interaction with transfer eps: creates secondaries and updates the track.
    virtual void produceSecondaries(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries, double transfer) = 0;

    // --------------------------------------------------------------------------
    // DiscreteProcess Overrides
    // --------------------------------------------------------------------------
    // 1 / (n * sigma(eps > eps_c)) from the lookup table.
    double getMeanFreePath(const Track& track, const Material& material) const override;

    // Optional discontinuity points in dsigma/deps (e.g., k_max' in bremsstrahlung)
    // where numerical integration is split into separate sub-intervals.
    virtual std::vector<double> getIntegrationBreakpoints(const Particle&, double, const Material&) const {return {};}

public:
    // --------------------------------------------------------------------------
    // Construction & Configuration
    // --------------------------------------------------------------------------
    EnergyLossProcess(std::string name, double cut)
        : DiscreteProcess(std::move(name)), productionCut(cut) {}

    // Must be set before the run: lookup tables depend on the production cut.
    void setProductionCut(double cut);

    double getProductionCut() const {return productionCut;}

    // --------------------------------------------------------------------------
    // Stepping Interface (DiscreteProcess overrides)
    // --------------------------------------------------------------------------
    // Does not limit the step; saves restricted dE/dx(T_pre) for alongStepDoIt.
    double alongStepGetPhysicalInteractionLength(const Track& track, const Material& material) override;

    // Continuous part: applies restricted loss via trapezoid rule between T_pre and T_now.
    void alongStepDoIt(Track& track, const Material& material, double stepLength) override;

    // Discrete part: accepts via integral approach, samples eps >= eps_c, calls produceSecondaries.
    void postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;

    // --------------------------------------------------------------------------
    // Physics & Validation Interface
    // Public so validation can compare against PDG tables and integrate dsigma/deps directly.
    // --------------------------------------------------------------------------
    // Differential cross section dsigma/deps per atom (cm^2 / energy), defined by derived class.
    virtual double computeDifferentialCrossSectionPerAtom(const Particle& particle,
                                                          double kineticEnergy, double transfer, const Material& material) const = 0;

    // Integrals of dsigma/deps over [lo, hi] (clipped to allowed range), per atom.
    double computeCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const;

    double computeEnergyLossPerAtom(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const;

    // Unrestricted and restricted stopping power (energy / cm).
    double computeDEDX(const Particle& particle, double kineticEnergy, const Material& material) const;

    double computeRestrictedDEDX(const Particle& particle, double kineticEnergy, const Material& material, double cut) const;

    // Allowed transfer range at kineticEnergy (public wrappers for validation).
    double getMinTransfer(const Particle& particle, const Material& material) const {
        return getMinEnergyTransfer(particle, material);
    }

    double getMaxTransfer(const Particle& particle, double kineticEnergy, const Material& material) const {
        return getMaxEnergyTransfer(particle, kineticEnergy, material);
    }

    // Samples eps from dsigma/deps on [cut, eps_max] directly, bypassing lookup tables.
    double sampleEnergyTransfer(const Particle& particle, double kineticEnergy, double cut, const Material& material) const;

    // --------------------------------------------------------------------------
    // Diagnostics
    // --------------------------------------------------------------------------
    long getMajorantViolations() const {return majorantViolations;}
    long getSigmaMaxViolations() const {return sigmaMaxViolations;}
};