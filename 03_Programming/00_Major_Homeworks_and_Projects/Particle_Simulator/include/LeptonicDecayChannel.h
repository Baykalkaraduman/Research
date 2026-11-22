#pragma once

#include "DecayChannel.h"

// ==============================================================================
// Class: LeptonicDecayChannel
// Three-body leptonic decay with the V-A matrix element of the Standard Model:
//     mu-  -> e-  anti_nu_e  nu_mu         mu+  -> e+  nu_e       anti_nu_mu
//     tau- -> l-  anti_nu_l  nu_tau        tau+ -> l+  nu_l       anti_nu_tau
//
// Daughter order: { charged lepton (1), neutrino of the lepton's flavour (2),
//                   neutrino of the parent's flavour (3) }
// Spin-averaged matrix element (Fermi theory, massless neutrinos):
//     |M|^2  ~  (P . p2) (p1 . p3)
// P = parent. In terms of the Dalitz variable s13 = (p1 + p3)^2 = (P - p2)^2:
//     P . p2  = (M^2 - s13) / 2,       p1 . p3 = (s13 - m1^2) / 2
// so |M|^2 depends on s13 only:   |M|^2 ~ (M^2 - s13)(s13 - m1^2).
//
// This gives the Michel spectrum of the charged lepton (x = E / E_max):
//     dN/dx ~ sqrt(x^2 - x0^2) (3x - 2x^2 - x0^2),   x0 = m1 / E_max
// which for m1 -> 0 is the familiar x^2 (3 - 2x). The mean is <x> = 0.7, so the
// electron takes 35 % of the muon's rest energy on average and the two
// neutrinos 65 %. (Before stage 4: "the electron gets 80 % of T", see Decay.cpp.)
// ==============================================================================
class LeptonicDecayChannel : public DecayChannel {
public:
    // daughters: { charged lepton, neutrino paired with the parent (2), neutrino paired with the lepton (3) }
    LeptonicDecayChannel(const Particle* parentDef, double ratio, const Particle* lepton,
                         const Particle* neutrino2, const Particle* neutrino3)
        : DecayChannel(parentDef, ratio, {lepton, neutrino2, neutrino3}) {}

    // -> LeptonicDecayChannel.cpp
    std::vector<LorentzVector> decayIt(double parentMass) const override;

    // The Michel spectrum above, not normalised (for the validation)
    static double michelSpectrum(double x, double x0) {return x > x0 ? std::sqrt(x * x - x0 * x0) * (3.0 * x - 2.0 * x * x - x0 * x0) : 0.0;}
};
