#include "GammaConversionMuon.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::alpha;
using PhysicalConstants::m_e;
using PhysicalConstants::m_mu;
using PhysicalConstants::r_e;

// ==============================================================================
// Constants & Helpers
// ==============================================================================
namespace {

constexpr double kSqrtE = 1.6487212707001282;

// Element-dependent parameters for the screening and nuclear form-factor formulas.
struct ElementParameters {
    double b = 0.0;
    double dn = 0.0;
    double zThird = 0.0;
    double wInfinity = 0.0;

    explicit ElementParameters(const Material& material) {
        const double z = material.getZ();
        const bool hydrogen = z < 1.5;
        b = hydrogen ? 202.4 : 183.0;
        dn = hydrogen ? 1.49 : 1.54 * std::pow(material.getA(), 0.27);
        zThird = 1.0 / std::cbrt(z); // Z^(-1/3)
        wInfinity = b * zThird * m_mu / (dn * m_e);
    }
};

} // namespace

// ==============================================================================
// Cross Section
// ==============================================================================
double GammaConversionToMuons::computeCrossSectionPerAtom(double photonEnergy, const Material& material) {
    if (photonEnergy <= 2.0 * m_mu || material.getZ() < 1.0) {
        return 0.0;
    }

    const ElementParameters el(material);
    const double z = material.getZ();
    const double rMu = r_e * m_e / m_mu;
    const double wMedium = 1.0 / (4.0 * el.dn * kSqrtE * m_mu);
    const double wSaturation = el.wInfinity / wMedium;
    const double powThreshold = 1.479 + 0.00799 * el.dn;
    const double powSaturation = -0.88;
    const double eCorrection = (-18.0 + 4347.0 / (el.b * el.zThird)) * Units::MeV;

    const double correction = 1.0 + 0.04 * std::log(1.0 + eCorrection / photonEnergy);
    const double eg = std::pow(1.0 - 4.0 * m_mu / photonEnergy, powThreshold)
                      - std::pow(std::pow(wSaturation, powSaturation) + std::pow(photonEnergy, powSaturation), 1.0 / powSaturation);

    return (7.0 / 9.0) * 4.0 * alpha * z * z * rMu * rMu * std::log(1.0 + wMedium * correction * eg);
}

// Below 5 m_mu: quadratic threshold scaling sigma(5 m_mu) * y^2
double GammaConversionToMuons::getCrossSectionPerAtom(double photonEnergy, const Material& material) {
    const double lowest = 2.0 * m_mu;
    const double limit = 5.0 * m_mu;
    if (photonEnergy <= lowest) {
        return 0.0;
    }
    if (photonEnergy < limit) {
        const double y = (photonEnergy - lowest) / (limit - lowest);
        return y * y * computeCrossSectionPerAtom(limit, material);
    }

    return computeCrossSectionPerAtom(photonEnergy, material);
}

double GammaConversionToMuons::getMeanFreePath(const Track& track, const Material& material) const {
    const double sigma = crossSectionFactor * material.getAtomDensity() * getCrossSectionPerAtom(track.getKineticEnergy(), material);
    return sigma > 0.0 ? 1.0 / sigma : GeoDefs::kInfinity;
}

// ==============================================================================
// Energy Sharing
// ==============================================================================
// Samples x+ uniformly in [x_min, 1 - x_min] (x_min = 1/2 - sqrt(1/4 - m/E)),
// accepted with:
//     (1 - 4/3 x+ x-) ln W(x+) / ln W_max,
//     W = W_inf (1 + (D_n sqrt(e) - 2) delta / m_mu) / (1 + B Z^(-1/3) sqrt(e) delta / m_e),
//     delta = m_mu^2 / (2 E x+ x-)    (minimum momentum transfer)
// Below 5 m_mu the muons share the energy equally (x+ = 1/2), as in Geant4.
double GammaConversionToMuons::sampleEnergyFraction(double photonEnergy, const Material& material) {
    const ElementParameters el(material);
    const double gammaMuonInv = m_mu / photonEnergy;
    const double xMin = (photonEnergy <= 5.0 * m_mu) ? 0.5 : 0.5 - std::sqrt(0.25 - gammaMuonInv);
    const double xMax = 1.0 - xMin;
    if (xMin >= 0.5) {
        return 0.5;
    }

    const double ds2 = el.dn * kSqrtE - 2.0;
    const double sBZ = kSqrtE * el.b * el.zThird / m_e;
    const double logWMaxInv = 1.0 / std::log(el.wInfinity * (1.0 + 2.0 * ds2 * gammaMuonInv) / (1.0 + 2.0 * sBZ * m_mu * gammaMuonInv));

    double xPlus = 0.0;
    double result = 0.0;
    do {
        xPlus = xMin + RandomEngine::uniNumGen() * (xMax - xMin);
        const double xPM = xPlus * (1.0 - xPlus);
        const double del = m_mu * m_mu / (2.0 * photonEnergy * xPM);
        const double w = el.wInfinity * (1.0 + ds2 * del / m_mu) / (1.0 + sBZ * del);
        const double xxp = 1.0 - (4.0 / 3.0) * xPM;
        result = (xxp > 0.0) ? xxp * std::log(w) * logWMaxInv : 0.0;
    } while (RandomEngine::uniNumGen() > result);

    return xPlus;
}

// ==============================================================================
// Secondary Production (photon -> mu+ mu-)
// ==============================================================================
// Samples Kelner angular variables t (polar, via u = sqrt(1/t - 1)), psi
// (momentum-transfer azimuth), and rho (momentum-transfer magnitude) via
// rejection, then constructs:
//     theta+ = (m/E)(u + rho cos(psi)/2) / x+,   theta- = (m/E)(u - rho cos(psi)/2) / x-
// with the two muons emitted at opposite azimuths phi0 +- phiHalf.
void GammaConversionToMuons::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    clearNumberOfInteractionLengthLeft();

    const double photonEnergy = track.getKineticEnergy();
    if (photonEnergy <= 2.0 * m_mu) {
        return;
    }

    const ElementParameters el(material);
    const double xPlus = sampleEnergyFraction(photonEnergy, material);
    const double xMinus = 1.0 - xPlus;
    const double xPM = xPlus * xMinus;
    const double gammaMuonInv = m_mu / photonEnergy;

    const double a027 = std::pow(material.getA(), 0.27);
    const double c1Num = 0.138 * a027;
    const double c2Term2 = m_e / (183.0 * el.zThird * m_mu);
    const double a3 = gammaMuonInv / (2.0 * xPM);
    const double a33 = a3 * a3;
    const double b1 = 1.0 / (4.0 * c1Num * c1Num);
    const double b3 = b1 * b1 * b1;
    const double a21 = a33 + b1;
    const double f1Max = -(1.0 - xPM) * (2.0 * b1 + (a21 + a33) * std::log(a33 / a21)) / (2.0 * b3);

    const int nMax = 1000;
    int n = 0;
    double thetaPlus = 0.0;
    double thetaMinus = 0.0;
    double phiHalf = 0.0;

    do {
        double t = 0.0;
        double f1 = 0.0;
        do {
            ++n;
            t = RandomEngine::uniNumGen();
            const double a34 = a33 / (t * t);
            const double a22 = a34 + b1;
            if (std::abs(b1) < 0.0001 * a34) {
                f1 = (1.0 - 2.0 * xPM + 4.0 * xPM * t * (1.0 - t)) / (12.0 * a34 * a34 * a34 * a34);
            } else {
                f1 = -(1.0 - 2.0 * xPM + 4.0 * xPM * t * (1.0 - t)) * (2.0 * b1 + (a22 + a34) * std::log(a34 / a22)) / (2.0 * b3);
            }
            if (f1 < 0.0 || f1 > f1Max) {
                f1 = 0.0; // should never happen (Geant4 prints a warning)
            }
        } while (RandomEngine::uniNumGen() * f1Max > f1 && n < nMax);

        const double f2Max = 1.0 - 2.0 * xPM * (1.0 - 4.0 * t * (1.0 - t));
        double psi = 0.0;
        double f2 = 0.0;
        do {
            ++n;
            psi = 2.0 * M_PI * RandomEngine::uniNumGen();
            f2 = 1.0 - 2.0 * xPM + 4.0 * xPM * t * (1.0 - t) * (1.0 + std::cos(2.0 * psi));
            if (f2 < 0.0 || f2 > f2Max) {
                f2 = 0.0;
            }
        } while (RandomEngine::uniNumGen() * f2Max > f2 && n < nMax);

        const double c2Term1 = gammaMuonInv / (2.0 * xPM * t);
        const double c22 = c2Term1 * c2Term1 + c2Term2 * c2Term2;
        const double c2 = 4.0 * c22 * c22 / std::sqrt(xPM);
        const double rhoMax = (1.0 / t - 1.0) * 1.9 / a027;
        const double beta = std::log((c2 + std::pow(rhoMax, 4)) / c2);
        const double rho = std::pow(c2 * (std::exp(beta * RandomEngine::uniNumGen()) - 1.0), 0.25);

        const double u = std::sqrt(1.0 / t - 1.0);
        const double xiHalf = 0.5 * rho * std::cos(psi);
        phiHalf = 0.5 * rho / u * std::sin(psi);
        thetaPlus = gammaMuonInv * (u + xiHalf) / xPlus;
        thetaMinus = gammaMuonInv * (u - xiHalf) / xMinus;

        if (n > nMax) {
            if (std::abs(thetaPlus) > M_PI) {
                thetaPlus = 0.0;
            }
            if (std::abs(thetaMinus) > M_PI) {
                thetaMinus = 0.0;
            }
        }
    } while (std::abs(thetaPlus) > M_PI || std::abs(thetaMinus) > M_PI);

    const double phi0 = 2.0 * M_PI * RandomEngine::uniNumGen();
    Vector3D plusDirection(std::sin(thetaPlus) * std::cos(phi0 + phiHalf), std::sin(thetaPlus) * std::sin(phi0 + phiHalf), std::cos(thetaPlus));
    Vector3D minusDirection(-std::sin(thetaMinus) * std::cos(phi0 - phiHalf), -std::sin(thetaMinus) * std::sin(phi0 - phiHalf), std::cos(thetaMinus));
    plusDirection.rotateUz(track.getDirection());
    minusDirection.rotateUz(track.getDirection());

    secondaries.emplace_back(muonPlusDef, track.getPosition(), plusDirection, xPlus * photonEnergy - m_mu);
    secondaries.emplace_back(muonMinusDef, track.getPosition(), minusDirection, xMinus * photonEnergy - m_mu);
    track.setKineticEnergy(0.0);
    track.killTrack();

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] mu+ mu- pair, x+ = " << xPlus << " (E_gamma = " << photonEnergy / Units::GeV << " GeV)\n";
    }
}