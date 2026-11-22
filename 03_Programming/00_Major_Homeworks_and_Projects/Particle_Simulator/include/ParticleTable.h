#pragma once

#include "DecayTable.h"
#include "Particle.h"

#include <memory>
#include <string>
#include <vector>

// ==============================================================================
// Class: ParticleTable
// Owns every particle definition of the simulation and their decay tables.
// Equivalent to G4ParticleTable (+ the G4MuonMinus, G4TauMinus, G4PionPlus, ...
// definition classes, which build their decay tables).
//
// Why a table now: a decay produces up to 6 different particle types (tau- ->
// pi- pi- pi+ pi0 nu_tau), and the channels need pointers to all of them. Before
// stage 4 every program created its own e-, gamma, ...; now there is one place
// that knows the particles and how each one decays.
//
// Contents:
//   e-  e+  gamma
//   mu-  mu+           -> e nu nu (100 %, V-A, Michel spectrum)
//   tau- tau+          -> e nu nu, mu nu nu, hadronic (see ParticleTable.cpp)
//   pi+  pi-           -> mu nu (99.988 %), e nu (0.012 %)
//   pi0                -> gamma gamma
//   nu_e anti_nu_e nu_mu anti_nu_mu nu_tau anti_nu_tau   (stable, no interactions)
// ==============================================================================
class ParticleTable {
private:
    std::vector<std::unique_ptr<Particle>> particles;
    std::vector<std::unique_ptr<DecayTable>> decayTables;

    // Create and store one particle (-> ParticleTable.cpp)
    Particle* insert(int pdg, const std::string& name, double mass, double charge, double spin, double lifetime);

    // Build the decay tables, after all particles exist (-> ParticleTable.cpp)
    void buildDecayTables();

public:
    // Creates all particles and decay tables (-> ParticleTable.cpp)
    ParticleTable();

    // The table hands out pointers to its particles: it must not be copied
    ParticleTable(const ParticleTable&) = delete;
    ParticleTable& operator=(const ParticleTable&) = delete;

    // nullptr if unknown (G4ParticleTable::FindParticle). -> ParticleTable.cpp
    const Particle* findParticle(const std::string& name) const;
    const Particle* findParticle(int pdg) const;

    int entries() const { return static_cast<int>(particles.size()); }
    const Particle* getParticle(int index) const { return particles[index].get(); }
};
