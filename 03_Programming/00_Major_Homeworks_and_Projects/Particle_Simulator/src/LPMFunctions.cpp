#include "LPMFunctions.h"
#include "PhysicalConstants.h"
#include "Units.h"

#include <algorithm>
#include <cmath>

using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

namespace {

constexpr double kAlphaEm = 1.0 / 137.035999084;
constexpr double kHbarCGeVCm = 1.973269804e-14; // hbar*c [GeV * cm]
constexpr double kLpmConstGeVPerCm = (m_e * m_e * kAlphaEm) / (4.0 * M_PI * kHbarCGeVCm); // ~7684.4 GeV/cm

} // namespace

// ==============================================================================
// Characteristic Energy Scales
// ==============================================================================
double LPMFunctions::computeLPMEnergy(const Material& material) {
    const double x0 = material.getRadLength(); // [cm]
    return (x0 > 0.0) ? kLpmConstGeVPerCm * x0 : 0.0;
}

double LPMFunctions::computePlasmaEnergy(const Material& material) {
    const double electronDensity = material.getAtomDensity() * material.getZ(); // [1/cm^3]
    if (electronDensity <= 0.0) {
        return 0.0;
    }
    return kHbarCGeVCm * std::sqrt(4.0 * M_PI * r_e * electronDensity);
}

// ==============================================================================
// Migdal Functions: s1, xi(s), s(s'), phi(s), G(s)
// ==============================================================================
double LPMFunctions::computeS1(int z) {
    if (z < 1) {
        return 1.0;
    }
    const double ratio = std::cbrt(static_cast<double>(z)) / 184.15;
    return ratio * ratio;
}

double LPMFunctions::computeXi(double s, double s1) {
    if (s >= 1.0 || s1 >= 1.0) {
        return 1.0;
    }
    if (s <= s1) {
        return 2.0;
    }
    return 1.0 + std::log(s) / std::log(s1);
}

// Solves the implicit equation s = s' / sqrt(xi(s)) via contractive fixed-point iteration
double LPMFunctions::solveMigdalS(double sPrime, double s1) {
    if (sPrime <= 0.0) {
        return 0.0;
    }
    if (sPrime >= 1.0) {
        return sPrime;
    }
    double s = sPrime;
    for (int i = 0; i < 8; ++i) {
        s = sPrime / std::sqrt(computeXi(s, s1));
    }
    return s;
}

// Exact Migdal series: phi(s) = 6s - 6*pi*s^2 + 24*s^3 * sum_{k=1}^inf 1 / ((k+s)^2 + s^2)
// Asymptote (s >= 2):  phi(s) = 1 - 1 / (84 * s^4)
double LPMFunctions::computePhi(double s) {
    if (s <= 0.0) {
        return 0.0;
    }
    const double s2 = s * s;
    if (s >= 2.0) {
        return 1.0 - 1.0 / (84.0 * s2 * s2);
    }

    constexpr int kTerms = 64;
    double sum = 0.0;
    for (int k = 1; k <= kTerms; ++k) {
        const double u = static_cast<double>(k) + s;
        sum += 1.0 / (u * u + s2);
    }
    sum += std::atan(s / (static_cast<double>(kTerms) + 0.5 + s)) / s;

    const double val = 6.0 * s - 6.0 * M_PI * s2 + 24.0 * s2 * s * sum;
    return std::clamp(val, 0.0, 1.0);
}

// Exact Migdal series: G(s) = 12*pi*s^2 - 48*s^3 * sum_{k=0}^inf 1 / ((k+0.5+s)^2 + s^2)
// Asymptote (s >= 2):  G(s) = 1 - 31 / (1344 * s^4)
double LPMFunctions::computeG(double s) {
    if (s <= 0.0) {
        return 0.0;
    }
    const double s2 = s * s;
    if (s >= 2.0) {
        return 1.0 - (31.0 / 1344.0) / (s2 * s2);
    }

    constexpr int kTerms = 64;
    double sum = 0.0;
    for (int k = 0; k < kTerms; ++k) {
        const double u = static_cast<double>(k) + 0.5 + s;
        sum += 1.0 / (u * u + s2);
    }
    sum += std::atan(s / (static_cast<double>(kTerms) + s)) / s;

    const double val = 12.0 * M_PI * s2 - 48.0 * s2 * s * sum;
    return std::clamp(val, 0.0, 1.0);
}

// ==============================================================================
// Suppression Ratios Relative to Complete-Screening Bethe-Heitler (in [0, 1])
// ==============================================================================
double LPMFunctions::computeDielectricFactor(double electronTotalEnergy, double photonEnergy,
                                             const Material& material) {
    if (electronTotalEnergy <= m_e || photonEnergy <= 0.0 || photonEnergy >= electronTotalEnergy) {
        return 1.0;
    }
    const double kp = (electronTotalEnergy / m_e) * computePlasmaEnergy(material);
    if (kp <= 0.0) {
        return 1.0;
    }

    const double ratio = kp / photonEnergy;
    const double gammaDielectric = 1.0 + ratio * ratio;
    if (gammaDielectric < 1.0 + 1.0e-5) {
        return 1.0;
    }

    const double y = photonEnergy / electronTotalEnergy;
    const double oneMinusY = 1.0 - y;
    const double twoTerm = 2.0 * (1.0 + oneMinusY * oneMinusY);
    const double bhDenom = y * y + twoTerm;
    const double deNum = (y * y) / (gammaDielectric * gammaDielectric) + twoTerm / gammaDielectric;
    return std::clamp(deNum / bhDenom, 0.0, 1.0);
}

double LPMFunctions::computeBremsstrahlungFactor(double electronTotalEnergy, double photonEnergy,
                                                 const Material& material, bool includeDielectric) {
    if (electronTotalEnergy <= m_e || photonEnergy <= 0.0 || photonEnergy >= electronTotalEnergy) {
        return 1.0;
    }
    const double eLpm = computeLPMEnergy(material);
    if (eLpm <= 0.0) {
        return includeDielectric ? computeDielectricFactor(electronTotalEnergy, photonEnergy, material) : 1.0;
    }

    const double y = photonEnergy / electronTotalEnergy;
    const double sPrime = std::sqrt((eLpm * photonEnergy) /
                                    (8.0 * electronTotalEnergy * (electronTotalEnergy - photonEnergy)));

    double gammaDielectric = 1.0;
    if (includeDielectric) {
        const double kp = (electronTotalEnergy / m_e) * computePlasmaEnergy(material);
        const double ratio = kp / photonEnergy;
        gammaDielectric = 1.0 + ratio * ratio;
    }

    if (sPrime > 5.0) {
        return includeDielectric ? computeDielectricFactor(electronTotalEnergy, photonEnergy, material) : 1.0;
    }

    const int z = static_cast<int>(std::lround(material.getZ()));
    const double s1 = computeS1(z);
    const double s = solveMigdalS(sPrime, s1);
    const double sEff = gammaDielectric * s;
    double xi = computeXi(sEff, s1);

    // Match Geant4: suppress Migdal's xi(s)*phi(s) > 1 overshoot
    const double phiS = computePhi(s);
    if (phiS > 0.0 && (xi * phiS > 1.0 || sPrime > 0.57)) {
        xi = 1.0 / phiS;
    }

    const double phi = (gammaDielectric == 1.0) ? phiS : computePhi(sEff);
    const double g = computeG(sEff);

    const double oneMinusY = 1.0 - y;
    const double twoTerm = 2.0 * (1.0 + oneMinusY * oneMinusY);
    const double bhDenom = y * y + twoTerm;
    const double lpmNum = xi * ((y * y / (gammaDielectric * gammaDielectric)) * g + (twoTerm / gammaDielectric) * phi);

    return std::clamp(lpmNum / bhDenom, 0.0, 1.0);
}

double LPMFunctions::computePairProductionFactor(double photonEnergy, double pairFraction,
                                                 const Material& material) {
    if (photonEnergy <= 2.0 * m_e || pairFraction <= 0.0 || pairFraction >= 1.0) {
        return 1.0;
    }
    const double eLpm = computeLPMEnergy(material);
    if (eLpm <= 0.0) {
        return 1.0;
    }

    const double ePlus = pairFraction * photonEnergy;
    const double eMinus = photonEnergy - ePlus;
    const double sPrime = std::sqrt((eLpm * photonEnergy) / (8.0 * ePlus * eMinus));
    if (sPrime > 5.0) {
        return 1.0;
    }

    const int z = static_cast<int>(std::lround(material.getZ()));
    const double s1 = computeS1(z);
    const double s = solveMigdalS(sPrime, s1);
    double xi = computeXi(s, s1);

    const double phi = computePhi(s);
    if (phi > 0.0 && (xi * phi > 1.0 || sPrime > 0.57)) {
        xi = 1.0 / phi;
    }

    const double twoTerm = 2.0 * (pairFraction * pairFraction + (1.0 - pairFraction) * (1.0 - pairFraction));
    const double bhDenom = 1.0 + twoTerm;
    const double lpmNum = xi * (computeG(s) + twoTerm * phi);

    return std::clamp(lpmNum / bhDenom, 0.0, 1.0);
}