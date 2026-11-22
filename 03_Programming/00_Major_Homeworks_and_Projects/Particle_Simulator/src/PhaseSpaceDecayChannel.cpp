#include "PhaseSpaceDecayChannel.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// ==============================================================================
// N-body phase space
// ==============================================================================
// Idea: build the decay as a chain of two-body decays
//     M = M_{N-1} -> M_{N-2} + m_{N-1},   M_{N-2} -> M_{N-3} + m_{N-2},  ...,  M_1 -> m_0 + m_1
// where M_k is the invariant mass of the first k+1 daughters. The M_k are chosen
// by dividing the kinetic energy T = M - sum(m_i) at N-2 uniform random points:
//     M_k = m_0 + ... + m_k + r_k T,    0 = r_0 < r_1 < ... < r_{N-1} = 1
// For uniform r's, the phase-space density of the configuration is proportional to
//     w = p*_1 p*_2 ... p*_{N-1},       p*_k = twoBodyMomentum(M_k, M_{k-1}, m_k)
// so we accept it with probability w / w_max (accept-reject). w_max is an upper
// bound: each p*_k is largest when M_k is largest and M_{k-1} smallest.
// The two-body case has one p* only, w = w_max, accepted at the first try.
//
// Then the momenta: start with daughters 0 and 1 back to back in the rest frame of
// M_1; for k = 2, 3, ...: in the rest frame of M_k, the subsystem (0..k-1) moves
// with momentum p*_k in a random direction and daughter k goes the opposite way.
// Boost the subsystem, add daughter k, repeat. After the last step we are in the
// parent's rest frame.
// ==============================================================================
std::vector<LorentzVector> PhaseSpaceDecayChannel::decayIt(double parentMass) const {
    const int n = getNumberOfDaughters();
    if (n < 2) {
        throw std::runtime_error("PhaseSpaceDecayChannel: needs at least two daughters (" + getDescription() + ")");
    }
    std::vector<double> mass(n);
    for (int i = 0; i < n; ++i) {
        mass[i] = daughters[i]->getMass();
    }
    const double kinetic = parentMass - getSumOfDaughterMasses();
    if (kinetic < 0.0) {
        throw std::runtime_error("PhaseSpaceDecayChannel: below threshold (" + getDescription() + ")");
    }

    // Upper bound of the weight
    double weightMax = 1.0;
    double massBelowMin = 0.0;       // smallest possible M_{k-1}
    double massAboveMax = kinetic;   // largest possible M_k
    massAboveMax += mass[0];
    for (int k = 1; k < n; ++k) {
        massBelowMin += mass[k - 1];
        massAboveMax += mass[k];
        weightMax *= twoBodyMomentum(massAboveMax, massBelowMin, mass[k]);
    }

    // Sample the intermediate masses until accepted
    std::vector<double> invariantMass(n);    // M_k
    std::vector<double> momentum(n, 0.0);    // p*_k, k >= 1
    std::vector<double> r(n);
    while (true) {
        r[0] = 0.0;
        r[n - 1] = 1.0;
        for (int k = 1; k < n - 1; ++k) {
            r[k] = RandomEngine::uniNumGen();
        }
        std::sort(r.begin() + 1, r.end() - 1);

        double massSum = 0.0;
        for (int k = 0; k < n; ++k) {
            massSum += mass[k];
            invariantMass[k] = massSum + r[k] * kinetic;
        }
        invariantMass[n - 1] = parentMass; // exactly, not up to rounding

        double weight = 1.0;
        for (int k = 1; k < n; ++k) {
            momentum[k] = twoBodyMomentum(invariantMass[k], invariantMass[k - 1], mass[k]);
            weight *= momentum[k];
        }
        if (weight >= RandomEngine::uniNumGen() * weightMax) {
            break;
        }
    }

    // Build the momenta, innermost two-body decay first
    std::vector<LorentzVector> products(n);
    const Vector3D dir = randomDirection();
    products[0] = LorentzVector(dir * momentum[1], std::sqrt(momentum[1] * momentum[1] + mass[0] * mass[0]));
    products[1] = LorentzVector(dir * (-momentum[1]), std::sqrt(momentum[1] * momentum[1] + mass[1] * mass[1]));
    for (int k = 2; k < n; ++k) {
        const Vector3D nk = randomDirection();
        const double q = momentum[k];
        // Subsystem (0..k-1), mass M_{k-1}, moves with momentum +q nk
        const double subEnergy = std::sqrt(q * q + invariantMass[k - 1] * invariantMass[k - 1]);
        const Vector3D beta = nk * (q / subEnergy);
        const double gamma = subEnergy / invariantMass[k - 1];
        for (int i = 0; i < k; ++i) {
            products[i].boost(beta, gamma);
        }
        products[k] = LorentzVector(nk * (-q), std::sqrt(q * q + mass[k] * mass[k]));
    }
    return products;
}
