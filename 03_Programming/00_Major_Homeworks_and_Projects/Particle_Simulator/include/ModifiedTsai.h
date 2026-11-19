#pragma once

// ==============================================================================
// Class: ModifiedTsai
// Polar angle of an e- / e+ produced in pair production, and of bremsstrahlung photon, relative to the parent direction.
// The physics: an ultra-relativistic particle radiates / is created inside a cone
// of opening ~ 1/gamma = m/E around the parent direction. Tsai's angular
// distribution is approximated by two terms in u = E theta / m:
//     f(u) ~ u exp(-a u) + 27 u exp(-3 a u),   a = 0.625
// The two terms have the weights 1/4 and 3/4 (their integrals are 1/a^2 and
// 3/a^2); each u exp(-c u) is sampled as the sum of two exponential numbers.
// ==============================================================================
class ModifiedTsai {
public:
    // cos(theta) for a particle of this kinetic energy (an electron's mass is used,
    // also for the photon of bremsstrahlung: the cone is set by the electron).
    static double sampleCosTheta(double kineticEnergy);
};
