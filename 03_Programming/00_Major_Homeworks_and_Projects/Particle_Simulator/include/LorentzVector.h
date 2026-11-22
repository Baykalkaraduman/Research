#pragma once

#include "Vector3D.h"

#include <cmath>

// ==============================================================================
// Struct: LorentzVector
// A four-momentum (p, E).
// ==============================================================================
struct LorentzVector {
    Vector3D p;        // momentum (GeV)
    double energy;     // total energy (GeV)

    LorentzVector(const Vector3D& momentum = Vector3D(), double totalEnergy = 0.0) : p(momentum), energy(totalEnergy) {}
    // ========================================================
    // Getters
    // ========================================================
    Vector3D vect() const {return p;}
    double E() const {return energy;}

    // Invariant mass squared m^2 = E^2 - p^2 (can be slightly negative from rounding)
    double m2() const {return energy * energy - p.dot(p);}
    double m() const {const double mm = m2(); return mm > 0.0 ? std::sqrt(mm) : 0.0;}

    // Velocity of the frame in which this four-vector is at rest: beta = p / E

    Vector3D boostVector() const {return p / energy;}

    LorentzVector operator+(const LorentzVector& v) const {return LorentzVector(p + v.p, energy + v.energy);}

    // Lorentz boost with velocity beta (|beta| < 1). gamma is passed separately
    // because 1/sqrt(1 - beta^2) loses digits for a fast parent: at beta*gamma = 1e4,
    // 1 - beta^2 = 1e-8 and gamma would keep only ~8 digits. The caller knows gamma = E/m exactly.
    LorentzVector& boost(const Vector3D& beta, double gamma);
};
