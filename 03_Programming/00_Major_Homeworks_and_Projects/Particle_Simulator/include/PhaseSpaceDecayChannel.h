#pragma once

#include "DecayChannel.h"

// ==============================================================================
// Class: PhaseSpaceDecayChannel
// Daughters distributed according to pure Lorentz-invariant phase space: every
// kinematically allowed configuration is equally likely (|M|^2 = constant).
//
//   2 bodies:  exact. Back to back, |p| = twoBodyMomentum(M, m1, m2), isotropic.
//              pi -> mu nu, pi0 -> gamma gamma, tau -> pi nu.
//   N bodies:  GENBOD (F. James, CERN 68-15), the algorithm of ROOT's TGenPhaseSpace
//
// ==============================================================================
class PhaseSpaceDecayChannel : public DecayChannel {
public:
    PhaseSpaceDecayChannel(const Particle* parentDef, double ratio, std::vector<const Particle*> daughterDefs)
        : DecayChannel(parentDef, ratio, std::move(daughterDefs)) {}

    // -> PhaseSpaceDecayChannel.cpp
    std::vector<LorentzVector> decayIt(double parentMass) const override;
};
