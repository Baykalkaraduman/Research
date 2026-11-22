#include "DecayChannel.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>

// ==============================================================================
// Helpers for all decay channels
// ==============================================================================

// cos(theta) uniform in [-1, 1], phi uniform in [0, 2 pi): isotropic
Vector3D DecayChannel::randomDirection() {
    const double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    return Vector3D(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

double DecayChannel::getSumOfDaughterMasses() const {
    double sum = 0.0;
    for (const Particle* d : daughters) {
        sum += d->getMass();
    }
    return sum;
}

std::string DecayChannel::getDescription() const {
    std::string text = parent->getName() + " ->";
    for (const Particle* d : daughters) {
        text += " " + d->getName();
    }
    return text;
}

// p* = sqrt[(M^2 - (m1 + m2)^2) (M^2 - (m1 - m2)^2)] / (2 M)
// Written as products of sums and differences: no M^4 - ... cancellation.
double DecayChannel::twoBodyMomentum(double parentMass, double m1, double m2) {
    const double a = (parentMass - m1 - m2) * (parentMass + m1 + m2);
    const double b = (parentMass - m1 + m2) * (parentMass + m1 - m2);
    if (a <= 0.0 || b <= 0.0) {
        return 0.0; // at (or below) threshold
    }
    return std::sqrt(a * b) / (2.0 * parentMass);
}
