#pragma once

// ==============================================================================
// Class: ModifiedMephi
// Polar angle of the photon of muon bremsstrahlung and of the e- / e+ of muon pair production, relative to the muon direction.
//
// The secondaries go into a cone of ~ 1/gamma of the MUON (gamma = E / M), not of
// the electron: the muon is the particle that radiates. With r = gamma * theta,
// the distribution is
//     f(r) dr ~ r dr / (1 + r^2)^2
// cut at r_max = gamma * pi/2 * min(1, gamma M / k - 1), k = energy of the
// secondary. Its cumulative is r^2 / (1 + r^2), so with x = u r_max^2 / (1 + r_max^2):
//     r = sqrt(x / (1 - x)),   theta = r / gamma
// One random number, no rejection.
// ==============================================================================
class ModifiedMephi {
public:
    // cos(theta) of a secondary of energy "secondaryEnergy" emitted by a particle
    // of mass `mass` and kinetic energy "kineticEnergy".
    static double sampleCosTheta(double kineticEnergy, double secondaryEnergy, double mass);
};
