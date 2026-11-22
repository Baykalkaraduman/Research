#pragma once

#include "LorentzVector.h"
#include "Particle.h"

#include <string>
#include <vector>

// ==============================================================================
// Abstract Class: DecayChannel
// Lists how particles can decay into daughters and their branching ratios
//
// decayIt() returns the daughters' four-momenta in the PARENT'S REST FRAME, in
// the same order as getDaughters(). The Decay process picks a channel, calls
// decayIt() and boosts the result into the lab.
//
// Derived classes only differ in HOW the momenta are distributed:
//   PhaseSpaceDecayChannel  pure phase space (two-body exact, N-body GENBOD)
//   LeptonicDecayChannel    V-A matrix element: mu -> e nu nu, tau -> l nu nu
// ==============================================================================
class DecayChannel {
protected:
    const Particle* parent;
    double branchingRatio;
    std::vector<const Particle*> daughters;

    // Isotropic unit vector.
    static Vector3D randomDirection();

public:
    DecayChannel(const Particle* parentDef, double ratio, std::vector<const Particle*> daughterDefs)
    : parent(parentDef), branchingRatio(ratio), daughters(std::move(daughterDefs)) {}
    virtual ~DecayChannel() = default;

    const Particle* getParent() const {return parent;}
    double getBR() const {return branchingRatio;}
    int getNumberOfDaughters() const {return static_cast<int>(daughters.size());}
    const Particle* getDaughter(int i) const {return daughters[i];}

    // Sum of the daughter masses, and the energy the decay releases
    // Q = M - sum(m_i): the kinetic energy the daughters share in the rest frame.
    double getSumOfDaughterMasses() const;
    double getQValue() const {return parent->getMass() - getSumOfDaughterMasses();}

    // "tau- -> pi- pi0 nu_tau" (for printing). -> DecayChannel.cpp
    std::string getDescription() const;

    // Momentum of each daughter in the rest frame of a two-body decay M -> m1 m2
    // (G4PhaseSpaceDecayChannel::Pmx, ROOT's TGenPhaseSpace PDK). -> DecayChannel.cpp
    static double twoBodyMomentum(double parentMass, double m1, double m2);

    // The daughters' four-momenta in the parent's rest frame
    virtual std::vector<LorentzVector> decayIt(double parentMass) const = 0;
};
