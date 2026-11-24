#pragma once

#include "PhysicsProcess.h"

// ==============================================================================
// Derived Class: MuonCapture
// Owns the stopped negative muon lifecycle in matter (Z >= 1), matching
// Geant4's G4MuonMinusCapture design:
//   1. Prompt muonic atomic cascade (K_alpha, L_alpha, local outer cascade = B_1s)
//   2. Combined lifetime: tau_eff = 1 / (Lambda_c(Z) + Q(Z) * Lambda_free)
//   3. Branching between nuclear capture ("muMinusCapture") and bound
//      decay-in-orbit ("Decay", with effective mass m_mu - B_1s).
// Capture rates and Z_eff: Suzuki, Measday & Roalsvig, Phys. Rev. C 35 (1987).
// ==============================================================================
class MuonCapture : public PhysicsProcess {
private:
    // --------------------------------------------------------------------------
    // Secondary Particle Definitions & State
    // --------------------------------------------------------------------------
    const Particle* gammaDef;
    const Particle* electronDef;
    const Particle* nuMuDef;     // nu_mu     (PDG  14)
    const Particle* antiNuEDef;  // anti_nu_e (PDG -12)
    const Particle* neutronDef;  // n         (PDG 2112, optional)
    bool emitMuonicXrays = true;

    // Bookkeeping for the most recent atRestDoIt call (for exact energy balance)
    double lastNuclearMassChange = 0.0; // DeltaM_nuc + sum(S_n) [GeV]
    int lastNeutronCount = 0;

    // --------------------------------------------------------------------------
    // Internal Helpers
    // --------------------------------------------------------------------------
    [[nodiscard]] static Vector3D sampleIsotropicDirection();
    [[nodiscard]] static double getFreeDecayRatePerNs(const Track& track);

    void emitAtomicCascade(int z, const Vector3D& pos, Step& step, std::vector<Track>& secondaries) const;
    void performBoundDecay(double binding1s, const Vector3D& pos, Step& step, std::vector<Track>& secondaries);
    void performNuclearCapture(int z, const Material& material, const Vector3D& pos, Step& step,
                               std::vector<Track>& secondaries);

public:
    // --------------------------------------------------------------------------
    // Construction & Configuration
    // --------------------------------------------------------------------------
    MuonCapture(const Particle* gamma, const Particle* electron, const Particle* nuMu,
                const Particle* antiNuE, const Particle* neutron = nullptr)
        : PhysicsProcess("muMinusCapture"), gammaDef(gamma), electronDef(electron),
          nuMuDef(nuMu), antiNuEDef(antiNuE), neutronDef(neutron) {}

    void setMuonicXrays(bool enable) {emitMuonicXrays = enable;}

    // --------------------------------------------------------------------------
    // Process Interface
    // --------------------------------------------------------------------------
    [[nodiscard]] bool isApplicable(const Particle& particle) const override {return particle.getPdgId() == 13;}
    [[nodiscard]] bool isAtRestApplicable(const Particle& particle, const Material& material) const override {
        return particle.getPdgId() == 13 && material.getZ() >= 1.0;
    }
    [[nodiscard]] bool isAtRestExclusive(const Particle& particle, const Material& material) const override {
        return isAtRestApplicable(particle, material);
    }

    [[nodiscard]] double atRestGetPhysicalInteractionLength(const Track& track, const Material& material) override;
    void atRestDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) override;

    // --------------------------------------------------------------------------
    // Last-Event Bookkeeping (for energy-balance & multiplicity checks)
    // --------------------------------------------------------------------------
    [[nodiscard]] double getLastNuclearMassChange() const {return lastNuclearMassChange;}
    [[nodiscard]] int getLastNeutronCount() const {return lastNeutronCount;}

    // --------------------------------------------------------------------------
    // Static Element Physics Lookups (public for validation)
    // --------------------------------------------------------------------------
    [[nodiscard]] static double getEffectiveZ(int z);
    [[nodiscard]] static double getHuffFactor(int z);
    [[nodiscard]] static double getCaptureRatePerNs(int z);
    [[nodiscard]] static double getCaptureProbability(int z);
    [[nodiscard]] static double getEffectiveLifetimeNs(int z);
    [[nodiscard]] static double get1sBindingEnergy(int z);
    [[nodiscard]] static double getKAlphaEnergy(int z);
    [[nodiscard]] static double getLAlphaEnergy(int z);
    [[nodiscard]] static double getNuclearMassChange(int z);
    [[nodiscard]] static double getFirstNeutronSeparationEnergy(int z);
};