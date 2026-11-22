#include "ParticleTable.h"
#include "LeptonicDecayChannel.h"
#include "PhaseSpaceDecayChannel.h"
#include "PhysicalConstants.h"

#include <stdexcept>

// ==============================================================================
// Particle Table 
// ==============================================================================

using namespace PhysicalConstants;

Particle* ParticleTable::insert(int pdg, const std::string& name, double mass, double charge, double spin, double lifetime) {
    particles.push_back(std::make_unique<Particle>(pdg, name, mass, charge, spin, lifetime));
    return particles.back().get();
}

ParticleTable::ParticleTable() {
    // PDG ID, name, mass, charge, spin, lifetime (-1 = stable)
    insert(11, "e-", m_e, -1.0, 0.5, -1.0);
    insert(-11, "e+", m_e, 1.0, 0.5, -1.0);
    insert(22, "gamma", 0.0, 0.0, 1.0, -1.0);

    insert(13, "mu-", m_mu, -1.0, 0.5, tau_mu);
    insert(-13, "mu+", m_mu, 1.0, 0.5, tau_mu);
    insert(15, "tau-", m_tau, -1.0, 0.5, tau_tau);
    insert(-15, "tau+", m_tau, 1.0, 0.5, tau_tau);

    insert(211, "pi+", m_pi, 1.0, 0.0, tau_pi);
    insert(-211, "pi-", m_pi, -1.0, 0.0, tau_pi);
    insert(111, "pi0", m_pi0, 0.0, 0.0, tau_pi0);

    // Neutrinos: massless, neutral, stable. No process is applicable to them
    // (Ionization and msc need a charge, the photon processes need PDG 22),
    // so they only move in straight lines until they leave the world. That is all
    // a neutrino does in a detector (interaction length in lead ~1e15 cm at 1 GeV).
    insert(12, "nu_e", 0.0, 0.0, 0.5, -1.0);
    insert(-12, "anti_nu_e", 0.0, 0.0, 0.5, -1.0);
    insert(14, "nu_mu", 0.0, 0.0, 0.5, -1.0);
    insert(-14, "anti_nu_mu", 0.0, 0.0, 0.5, -1.0);
    insert(16, "nu_tau", 0.0, 0.0, 0.5, -1.0);
    insert(-16, "anti_nu_tau", 0.0, 0.0, 0.5, -1.0);

    buildDecayTables();
}

const Particle* ParticleTable::findParticle(const std::string& name) const {
    for (const auto& p : particles) {
        if (p->getName() == name) return p.get();
    }
    return nullptr;
}

const Particle* ParticleTable::findParticle(int pdg) const {
    for (const auto& p : particles) {
        if (p->getPdgId() == pdg) return p.get();
    }
    return nullptr;
}

// ==============================================================================
// Decay tables
// ==============================================================================
// For each unstable particle X and its antiparticle the channels are the charge
// conjugates of each other: swap every particle with its antiparticle.
void ParticleTable::buildDecayTables() {
    auto get = [this](const std::string& name) {
        const Particle* p = findParticle(name);
        if (p == nullptr) throw std::runtime_error("ParticleTable: unknown particle " + name);
        return p;
    };
    auto newTable = [this](const std::string& name) {
        decayTables.push_back(std::make_unique<DecayTable>());
        // const_cast: the table was created above as non-const, findParticle only
        // hands out const pointers
        const_cast<Particle*>(findParticle(name))->setDecayTable(decayTables.back().get());
        return decayTables.back().get();
    };

    // ---------------------------------------------------------------------------
    // Muon: mu- -> e- anti_nu_e nu_mu  (100 %, G4MuonDecayChannel)
    // ---------------------------------------------------------------------------
    newTable("mu-")->insert(std::make_unique<LeptonicDecayChannel>(get("mu-"), 1.0, get("e-"), get("anti_nu_e"), get("nu_mu")));
    newTable("mu+")->insert(std::make_unique<LeptonicDecayChannel>(get("mu+"), 1.0, get("e+"), get("nu_e"), get("anti_nu_mu")));

    // ---------------------------------------------------------------------------
    // Charged pion: two-body, phase space is exact (spin 0)
    // ---------------------------------------------------------------------------
    // pi -> e nu is 1.23e-4, not 1/2: helicity suppression ~ (m_e / m_mu)^2.
    // Kept because it is a famous number and costs nothing.
    const double brPiE = 1.2327e-4;
    DecayTable* piPlus = newTable("pi+");
    piPlus->insert(std::make_unique<PhaseSpaceDecayChannel>(get("pi+"), 1.0 - brPiE, std::vector<const Particle*>{get("mu+"), get("nu_mu")}));
    piPlus->insert(std::make_unique<PhaseSpaceDecayChannel>(get("pi+"), brPiE, std::vector<const Particle*>{get("e+"), get("nu_e")}));
    DecayTable* piMinus = newTable("pi-");
    piMinus->insert(std::make_unique<PhaseSpaceDecayChannel>(get("pi-"), 1.0 - brPiE, std::vector<const Particle*>{get("mu-"), get("anti_nu_mu")}));
    piMinus->insert(std::make_unique<PhaseSpaceDecayChannel>(get("pi-"), brPiE, std::vector<const Particle*>{get("e-"), get("anti_nu_e")}));

    // Neutral pion: gamma gamma 98.8 %. The Dalitz decay e+ e- gamma (1.2 %) is
    // left out (it is not phase space: it peaks at small e+e- mass), so 100 %.
    newTable("pi0")->insert(std::make_unique<PhaseSpaceDecayChannel>(get("pi0"), 1.0, std::vector<const Particle*>{get("gamma"), get("gamma")}));

    // ---------------------------------------------------------------------------
    // Tau (PDG 2024). Leptonic
    // Hadronic: option 2 of the plan: real pions, which then ionise, scatter and
    // decay (pi0 -> gamma gamma showers), but have no strong interaction.
    // Momenta from phase space (no rho / a1 resonances, see PhaseSpaceDecayChannel.h).
    // ---------------------------------------------------------------------------
    const double brE = 0.1782, brMu = 0.1739;
    struct Hadronic { double br; int nPiZero; bool threeProng; };
    // threeProng: pi- pi- pi+ (instead of a single pi-)
    const std::vector<Hadronic> hadronic = {
        {0.1082, 0, false},  // pi- nu
        {0.2549, 1, false},  // pi- pi0 nu
        {0.0926, 2, false},  // pi- 2pi0 nu
        {0.0104, 3, false},  // pi- 3pi0 nu
        {0.0902, 0, true},   // pi- pi- pi+ nu
        {0.0449, 1, true},   // pi- pi- pi+ pi0 nu
    };
    // These six add up to 60.1 %. The rest of the hadronic width (4.7 %: kaon
    // channels, eta, ...) has no particles here; it is shared out among the pion
    // channels in proportion, so the hadronic total is 1 - BR(e) - BR(mu) = 64.8 %.
    double hadronicSum = 0.0;
    for (const Hadronic& h : hadronic) hadronicSum += h.br;
    const double scale = (1.0 - brE - brMu) / hadronicSum;

    for (int sign : {-1, +1}) {
        // sign = -1: tau-;  sign = +1: tau+ (all charges and particles conjugated)
        const bool minus = sign < 0;
        const Particle* tau = get(minus ? "tau-" : "tau+");
        const Particle* nuTau = get(minus ? "nu_tau" : "anti_nu_tau");
        const Particle* piSame = get(minus ? "pi-" : "pi+");      // same charge as the tau
        const Particle* piOpposite = get(minus ? "pi+" : "pi-");
        DecayTable* table = newTable(tau->getName());

        table->insert(std::make_unique<LeptonicDecayChannel>(tau, brE, get(minus ? "e-" : "e+"),
                                                             get(minus ? "anti_nu_e" : "nu_e"), nuTau));
        table->insert(std::make_unique<LeptonicDecayChannel>(tau, brMu, get(minus ? "mu-" : "mu+"),
                                                             get(minus ? "anti_nu_mu" : "nu_mu"), nuTau));
        for (const Hadronic& h : hadronic) {
            std::vector<const Particle*> d;
            if (h.threeProng) {
                d = {piSame, piSame, piOpposite};
            } else {
                d = {piSame};
            }
            for (int i = 0; i < h.nPiZero; ++i) d.push_back(get("pi0"));
            d.push_back(nuTau);
            table->insert(std::make_unique<PhaseSpaceDecayChannel>(tau, h.br * scale, d));
        }
    }
}
