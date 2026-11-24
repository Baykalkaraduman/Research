#include "MuPairProduction.h"
#include "ModifiedMephi.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using PhysicalConstants::alpha;
using PhysicalConstants::m_e;
using PhysicalConstants::r_e;

namespace {

const double sqrtE = 1.6487212707001282; // sqrt(e)

// Gauss-Legendre on [0, 1], 8 points
const double xgi[8] = {0.0198550717512320, 0.1016667612931865, 0.2372337950418355, 0.4082826787521750,
                       0.5917173212478250, 0.7627662049581645, 0.8983332387068135, 0.9801449282487680};
const double wgi[8] = {0.0506142681451880, 0.1111905172266870, 0.1568533229389435, 0.1813418916891810,
                       0.1813418916891810, 0.1568533229389435, 0.1111905172266870, 0.0506142681451880};

// Post-Born Coulomb correction spectrum F(y) = (1 - y) F_1(z) + y^2 F_2(z)
// with z = (M / m_e)^2 y^2 / (1 - y) (Ivanov, Schiller, Serbo 1998, Eqs. 20-23;
// Groom, Mokhov, Striganov 2001).
double computeCoulombSpectrum(double z, double y) {
    if (z <= 0.0) {
        return 0.0;
    }

    const double oneMinusY = 1.0 - y;
    const double y2 = y * y;
    const double lnZ = std::log(z);

    // Large-z asymptotic branch (Eq. 23)
    if (z > 1.0e4) {
        const double f1 = ((2.0 - 22.0 / (3.0 * z)) * lnZ + 6.0 + 5.0 / (9.0 * z)) / z;
        const double f2 = ((2.0 + 1.0 / z) * lnZ + 1.0 + 1.0 / (2.0 * z)) / z;
        return oneMinusY * f1 + y2 * f2;
    }

    // Small-z asymptotic branch (Eq. 22)
    if (z < 1.0e-3) {
        const double f = (7.0 / 3.0 + 8.0 * z / 15.0) * (-lnZ) + 20.0 / 9.0 + 511.0 * z / 450.0;
        return oneMinusY * f;
    }

    // Exact closed form (Eq. 21)
    const double sqrtZ = std::sqrt(z);
    const double sqrtFourPlusZ = std::sqrt(4.0 + z);
    const double s = (sqrtFourPlusZ / sqrtZ) * std::log(0.5 * (sqrtFourPlusZ + sqrtZ));

    const double f1 = 44.0 / (15.0 * z) - 16.0 / 15.0 - (7.0 / 3.0 + 8.0 * z / 15.0) * lnZ
                      + (-44.0 / (15.0 * z) + 4.0 / (4.0 + z) + 38.0 / 15.0 + 16.0 * z / 15.0) * s;

    const double f2 = -4.0 / (3.0 * z) - (7.0 / 6.0) * lnZ
                      + (-2.0 / (3.0 * z) + 8.0 / (z * (4.0 + z)) + 7.0 / 3.0) * s;

    return oneMinusY * f1 + y2 * f2;
}

} // namespace

double MuPairProduction::getMinEnergyTransfer(const Particle&, const Material&) const {
    return 4.0 * m_e;
}

// The residual energy E - eps must stay above (3/4) sqrt(e) Z^(1/3) M: below it the
// screening logarithms of the formula turn negative (G4MuPairProductionModel::
// MaxSecondaryEnergyForElement)
double MuPairProduction::getMaxEnergyTransfer(const Particle& particle, double kineticEnergy, const Material& material) const {
    return kineticEnergy + particle.getMass() * (1.0 - 0.75 * sqrtE * std::cbrt(material.getZ()));
}

// ==============================================================================
// dsigma/deps per atom
// ==============================================================================
// Variables of Kokoulin's formula:
//   beta = eps^2 / (2 E (E - eps))          (not a velocity)
//   xi   = (M/m_e)^2 beta (1 - rho^2) / 2   (the virtual photon's "mass" in units of M)
// and for the electron (e) and muon (m) diagrams
//   B_e, B_m   : the rho-dependent bracket functions (two regimes of xi each)
//   Y_e, Y_m   : effective screening parameters
//   Phi_e = B_e [ ln(B Z^(-1/3) sqrt((1+xi)(1+Y_e)) / (1 + s (1+xi)(1+Y_e)/(1-rho^2)))
//                 - 1/2 ln(1 + (3 m_e Z^(1/3) / 2M)^2 (1+xi)(1+Y_e)) ]
//   Phi_mu = B_m ln( (M/m_e) B Z^(-1/3) / ((3/2) Z^(1/3) (1 + s (1+xi)(1+Y_m)/(1-rho^2))) )
//   s = 2 sqrt(e) m_e B / (Z^(1/3) eps)    (screening)
// Z^(1/3) inside the logs is the nuclear-size correction.
// Includes the post-Born Coulomb correction of Ivanov, Schiller, Serbo (1998).
double MuPairProduction::computeDifferentialCrossSectionPerAtom(const Particle& particle, double kineticEnergy,
                                                                double transfer, const Material& material) const {
    if (transfer <= 4.0 * m_e) {
        return 0.0;
    }
    const double z = material.getZ();
    const double z13 = std::cbrt(z), z23 = z13 * z13;
    const double mass = particle.getMass();
    const double totalEnergy = kineticEnergy + mass;
    const double residualEnergy = totalEnergy - transfer;
    if (residualEnergy <= 0.75 * sqrtE * z13 * mass) {
        return 0.0;
    }

    // Integration variable t = ln(1 - rho), from tmn = ln(1 - rho_max) to 0
    const double a0 = 1.0 / (totalEnergy * residualEnergy);
    const double alf = 4.0 * m_e / transfer;
    const double rt = std::sqrt(1.0 - alf);
    const double delta = 6.0 * mass * mass * a0;
    const double tmnExp = alf / (1.0 + rt) + delta * rt; // = 1 - rho_max
    if (tmnExp >= 1.0) {
        return 0.0;
    }
    const double tmn = std::log(tmnExp);

    const double massRatio = mass / m_e;
    const double massRatio2 = massRatio * massRatio;

    // zeta: the atomic-electron term, Z (Z + zeta)
    const bool hydrogen = z < 1.5;
    const double bbb = hydrogen ? 202.4 : 183.0;
    const double g1 = hydrogen ? 4.4e-5 : 1.95e-5;
    const double g2 = hydrogen ? 4.8e-5 : 5.3e-5;
    double zeta = 0.0;
    const double z1Exp = totalEnergy / (mass + g1 * z23 * totalEnergy);
    if (z1Exp > 35.221047195922) { // i.e. 0.073 ln(z1Exp) - 0.26 > 0
        const double z2Exp = totalEnergy / (mass + g2 * z13 * totalEnergy);
        zeta = (0.073 * std::log(z1Exp) - 0.26) / (0.058 * std::log(z2Exp) - 0.14);
    }
    const double z2 = z * (z + zeta);
    const double screen0 = 2.0 * m_e * sqrtE * bbb / (z13 * transfer);
    const double beta = 0.5 * transfer * transfer * a0;
    const double xi0 = 0.5 * massRatio2 * beta;

    double sum = 0.0;
    for (int i = 0; i < 8; ++i) {
        const double rho = std::exp(tmn * xgi[i]) - 1.0; // rho in [-rho_max, 0]: the formula is even in rho
        const double rho2 = rho * rho;
        const double xi = xi0 * (1.0 - rho2);
        const double xi1 = 1.0 + xi;
        const double xii = 1.0 / xi;

        const double yeu = (4.0 * beta + 5.0) + (4.0 * beta - 1.0) * rho2;
        const double yed = (6.0 * beta + 2.0) * std::log(3.0 + xii) + (2.0 * beta - 1.0) * rho2 - 4.0 * beta;
        const double ymu = (6.0 * beta + 2.0) * (1.0 + rho2) + 6.0;
        const double ymd = (4.0 * beta + 3.0) * (1.0 + rho2) * std::log(3.0 + xi) + 2.0 - 3.0 * rho2;
        const double ye1 = 1.0 + yeu / yed;
        const double ym1 = 1.0 + ymu / ymd;

        double be, bm;
        if (xi <= 1000.0) {
            be = ((2.0 + rho2) * (1.0 + beta) + xi * (3.0 + rho2)) * std::log(1.0 + xii) + (1.0 - rho2 - beta) / xi1 - (3.0 + rho2);
        } else {
            be = 0.5 * (3.0 - rho2 + 2.0 * beta * (1.0 + rho2)) * xii; // large-xi expansion
        }
        if (xi >= 0.001) {
            const double a10 = (1.0 + 2.0 * beta) * (1.0 - rho2);
            bm = ((1.0 + rho2) * (1.0 + 1.5 * beta) - a10 * xii) * std::log(xi1) + xi * (1.0 - rho2 - beta) / xi1 + a10;
        } else {
            bm = 0.5 * (5.0 - rho2 + beta * (3.0 + rho2)) * xi;      // small-xi expansion
        }

        const double screen = screen0 * xi1 / (1.0 - rho2);
        const double ale = std::log(bbb / z13 * std::sqrt(xi1 * ye1) / (1.0 + screen * ye1));
        const double cre = 0.5 * std::log(1.0 + 2.25 * z23 * xi1 * ye1 / massRatio2);
        const double fe = std::max(0.0, (ale - cre) * be);
        const double fm = std::max(0.0, std::log(bbb * massRatio / (1.5 * z23 * (1.0 + screen * ym1))) * bm) / massRatio2;

        sum += wgi[i] * (1.0 + rho) * (fe + fm); // (1 + rho) = d rho / dt
    }

    const double factorForCross = alpha * alpha * r_e * r_e * 4.0 / (3.0 * M_PI);
    const double born = -tmn * sum * factorForCross * z2 * residualEnergy / (totalEnergy * transfer);

    // Post-Born Coulomb correction (Ivanov, Schiller, Serbo 1998, Eqs. 20-21)
    const double y = transfer / totalEnergy;
    const double zCoul = massRatio2 * transfer * transfer * a0; // = (M / m_e)^2 y^2 / (1 - y)
    const double fCoul = material.computeCoulombCorrection();
    const double coulomb = -factorForCross * z * z * fCoul * computeCoulombSpectrum(zCoul, y) / transfer;

    return std::max(0.0, born + coulomb);
}

// ==============================================================================
// e- and e+ (G4MuPairProductionModel::SampleSecondaries)
// ==============================================================================
//   r uniform in [-rho_max, rho_max],  E(e-) = (1 - r) eps / 2,  E(e+) = eps - E(e-)
// (total energies). Directions: ModifiedMephi with the pair energy, opposite
// azimuths. Muon: T - eps, direction from momentum conservation.
void MuPairProduction::produceSecondaries(Track& track, const Material&, Step&, std::vector<Track>& secondaries,
                                          double transfer) {
    const double kineticEnergy = track.getKineticEnergy();
    const double mass = track.getMass();
    const double totalEnergy = kineticEnergy + mass;

    const double rMax = (1.0 - 6.0 * mass * mass / (totalEnergy * (totalEnergy - transfer))) * std::sqrt(1.0 - 4.0 * m_e / transfer);
    const double r = rMax * (2.0 * RandomEngine::uniNumGen() - 1.0);
    const double electronTotal = 0.5 * (1.0 - r) * transfer;
    const double positronTotal = transfer - electronTotal;

    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    const double cosE = ModifiedMephi::sampleCosTheta(kineticEnergy, transfer, mass);
    const double sinE = std::sqrt((1.0 - cosE) * (1.0 + cosE));
    const double cosP = ModifiedMephi::sampleCosTheta(kineticEnergy, transfer, mass);
    const double sinP = std::sqrt((1.0 - cosP) * (1.0 + cosP));
    Vector3D electronDirection(sinE * std::cos(phi), sinE * std::sin(phi), cosE);
    Vector3D positronDirection(-sinP * std::cos(phi), -sinP * std::sin(phi), cosP);
    electronDirection.rotateUz(track.getDirection());
    positronDirection.rotateUz(track.getDirection());

    const double electronKinetic = std::max(electronTotal - m_e, 0.0);
    const double positronKinetic = std::max(positronTotal - m_e, 0.0);
    secondaries.emplace_back(electronDef, track.getPosition(), electronDirection, electronKinetic);
    secondaries.emplace_back(positronDef, track.getPosition(), positronDirection, positronKinetic);

    const double pe = std::sqrt(electronKinetic * (electronKinetic + 2.0 * m_e));
    const double pp = std::sqrt(positronKinetic * (positronKinetic + 2.0 * m_e));
    const Vector3D newMomentum = track.getDirection() * track.getMomentumMag() - electronDirection * pe - positronDirection * pp;
    track.deductEnergy(transfer);
    if (track.getKineticEnergy() > 0.0) {
        track.setDirection(newMomentum);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] e+e- pair of " << transfer / Units::GeV << " GeV, mu "
                  << kineticEnergy / Units::GeV << " -> " << track.getKineticEnergy() / Units::GeV << " GeV\n";
    }
}