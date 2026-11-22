#include "LorentzVector.h"

// ==============================================================================
// Lorentz boost
// ==============================================================================
// Split p into the part along beta and the part across it. Only the part along beta and the energy change:
//     p_par' = gamma (p_par + beta E)
//     E'     = gamma (E + beta . p)
// Written without splitting p:
//     p' = p + [ (gamma - 1) (beta . p) / beta^2  +  gamma E ] beta
LorentzVector& LorentzVector::boost(const Vector3D& beta, double gamma) {
    const double beta2 = beta.dot(beta);
    if (beta2 <= 0.0) {
        return *this; // no motion, nothing to do
    }
    const double betaDotP = beta.dot(p);
    // (gamma - 1) / beta^2 = gamma^2 / (gamma + 1): the right side does not
    // subtract two nearly equal numbers when beta is small
    const double gamma2 = gamma * gamma / (gamma + 1.0);
    p = p + beta * (gamma2 * betaDotP + gamma * energy);
    energy = gamma * (energy + betaDotP);
    return *this;
}
