#include "LeptonicDecayChannel.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// ==============================================================================
// Three-body decay: Dalitz plot + V-A weight
// ==============================================================================
// A three-body decay has 2 free variables once the orientation is ignored. The
// Dalitz variables s12 = (p1 + p2)^2 and s13 = (p1 + p3)^2 are the convenient
// ones: phase space is FLAT in (s12, s13). So:
//   1. draw (s12, s13) uniformly in the bounding box,
//   2. reject points outside the kinematic boundary,
//   3. accept with probability |M|^2 / |M|^2_max   (|M|^2 depends on s13 only),
//   4. turn (s12, s13) into energies and the angle between p1 and p2,
//   5. orient the whole event randomly in space (no polarisation).
//
// Energies in the rest frame:  E_i = (M^2 + m_i^2 - s_jk) / (2M)   ({i,j,k} = {1,2,3})
// Momentum conservation makes p1, p2, p3 a closed triangle, so
//     cos(theta_12) = (p3^2 - p1^2 - p2^2) / (2 p1 p2),
// and the point is inside the Dalitz plot exactly when |cos(theta_12)| <= 1.
std::vector<LorentzVector> LeptonicDecayChannel::decayIt(double parentMass) const {
    const double M = parentMass;
    const double m1 = daughters[0]->getMass();
    const double m2 = daughters[1]->getMass();
    const double m3 = daughters[2]->getMass();
    if (M < m1 + m2 + m3) {
        throw std::runtime_error("LeptonicDecayChannel: below threshold (" + getDescription() + ")");
    }
    const double M2 = M * M;

    const double s12Min = (m1 + m2) * (m1 + m2), s12Max = (M - m3) * (M - m3);
    const double s13Min = (m1 + m3) * (m1 + m3), s13Max = (M - m2) * (M - m2);
    const double sumOfSquares = M2 + m1 * m1 + m2 * m2 + m3 * m3;  // s12 + s13 + s23

    // |M|^2 ~ (M^2 + m2^2 - s13)(s13 - m1^2 - m3^2): a parabola in s13 with its
    // maximum at s13 = (M^2 + m2^2 + m1^2 + m3^2) / 2
    auto matrixElement = [&](double s13) { return (M2 + m2 * m2 - s13) * (s13 - m1 * m1 - m3 * m3); };
    const double halfWidth = 0.5 * (M2 + m2 * m2 - m1 * m1 - m3 * m3);
    const double matrixElementMax = halfWidth * halfWidth;

    double e1 = 0.0, e2 = 0.0, p1 = 0.0, p2 = 0.0, cos12 = 0.0;
    while (true) {
        const double s12 = s12Min + (s12Max - s12Min) * RandomEngine::uniNumGen();
        const double s13 = s13Min + (s13Max - s13Min) * RandomEngine::uniNumGen();
        const double s23 = sumOfSquares - s12 - s13;

        e1 = (M2 + m1 * m1 - s23) / (2.0 * M);
        e2 = (M2 + m2 * m2 - s13) / (2.0 * M);
        const double e3 = (M2 + m3 * m3 - s12) / (2.0 * M);
        if (e1 < m1 || e2 < m2 || e3 < m3) continue;

        p1 = std::sqrt(e1 * e1 - m1 * m1);
        p2 = std::sqrt(e2 * e2 - m2 * m2);
        const double p3sq = e3 * e3 - m3 * m3;
        if (p1 <= 0.0 || p2 <= 0.0) continue;
        cos12 = (p3sq - p1 * p1 - p2 * p2) / (2.0 * p1 * p2);
        if (std::abs(cos12) > 1.0) continue;  // outside the Dalitz plot

        if (matrixElement(s13) >= RandomEngine::uniNumGen() * matrixElementMax) break;
    }

    // Daughter 1 along a random direction, daughter 2 at theta_12 from it with a
    // random azimuth, daughter 3 balances the momentum
    const Vector3D n1 = randomDirection();
    const double sin12 = std::sqrt(std::max(0.0, 1.0 - cos12 * cos12));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    Vector3D n2(sin12 * std::cos(phi), sin12 * std::sin(phi), cos12);
    n2.rotateUz(n1);

    const Vector3D mom1 = n1 * p1;
    const Vector3D mom2 = n2 * p2;
    const Vector3D mom3 = (mom1 + mom2) * -1.0;
    // E3 from energy conservation (exact), not from |p3| (would differ by rounding)
    return {LorentzVector(mom1, e1), LorentzVector(mom2, e2), LorentzVector(mom3, M - e1 - e2)};
}
