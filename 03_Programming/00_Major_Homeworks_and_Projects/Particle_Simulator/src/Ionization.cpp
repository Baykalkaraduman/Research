#include "Ionization.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

// ROOT's Landau functions (CERNLIB DENLAN / DISLAN / RANLAN, in MathCore)
#include <Math/PdfFuncMathCore.h>    // ROOT::Math::landau_pdf
#include <Math/ProbFuncMathCore.h>   // ROOT::Math::landau_cdf
#include <Math/QuantFuncMathCore.h>  // ROOT::Math::landau_quantile

#include <algorithm> // std::max, std::min, std::upper_bound
#include <cmath>
#include <iostream>

using PhysicalConstants::K;
using PhysicalConstants::m_e;

// ==============================================================================
// Stopping power: model selection and low-energy extrapolation
// ==============================================================================
double Ionization::computeDEDX(const Particle& particle, double kineticEnergy, const Material& material) const {
    if (material.getDensity() <= 0.0 || material.getZ() <= 0.0) {
        return 0.0; // vacuum: no electrons to collide with
    }

    // Below the validity limit of the formulas, Bethe-Bloch's logarithm shrinks and
    // eventually goes NEGATIVE (for a muon in lead at beta*gamma < ~0.03). A negative
    // dE/dx would GIVE energy to the particle and it would never stop.
    // Instead: dE/dx(T) = dE/dx(T_low) * sqrt(T / T_low) below T_low. This is what
    // Geant4 does below the lowest energy of its tables, and it has the right shape
    // (the stopping power of a very slow ion falls like its velocity, ~ sqrt(T)).
    const double lowLimit = getLowEnergyLimit(particle);
    double energy = kineticEnergy;
    double scale = 1.0;
    if (energy < lowLimit) {
        scale = std::sqrt(energy / lowLimit);
        energy = lowLimit;
    }

    const double dedx = isElectronOrPositron(particle)
                      ? computeElectronDEDX(particle, energy, material)
                      : computeBetheBlochDEDX(particle, energy, material);
    return std::max(dedx, 0.0) * scale;
}

// Heavy particles: beta*gamma = 0.1, i.e. T = M (sqrt(1.01) - 1) ~ 0.005 M
// (0.53 MeV for a muon). Below this, shell corrections and the Barkas effect,
// which we do not have, change dE/dx by tens of percent.
// Electrons: 10 keV.
double Ionization::getLowEnergyLimit(const Particle& particle) const {
    if (isElectronOrPositron(particle)) {
        return 1.0 * Units::keV;
    }
    return particle.getMass() * (std::sqrt(1.0 + 0.1 * 0.1) - 1.0);
}

// ==============================================================================
// Maximum energy transfer to one electron in one collision
// ==============================================================================
//   heavy particle (PDG eq. 34.4):
//       Tmax = 2 m_e beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)
//   e- on e- (Moller): the two outgoing electrons are identical; by convention
//       the faster one is "the primary", so the other gets at most T/2.
//   e+ on e- (Bhabha): the particles are different, the positron can give all of T.
double Ionization::computeMaxSecondaryEnergy(const Particle& particle, double kineticEnergy) const {
    if (particle.getPdgId() == 11) {
        return 0.5 * kineticEnergy;
    }
    if (particle.getPdgId() == -11) {
        return kineticEnergy;
    }
    const double mass = particle.getMass();
    const double gamma = 1.0 + kineticEnergy / mass;
    const double betaGamma2 = kineticEnergy * (kineticEnergy + 2.0 * mass) / (mass * mass);
    const double ratio = m_e / mass;
    return 2.0 * m_e * betaGamma2 / (1.0 + 2.0 * gamma * ratio + ratio * ratio);
}

// ==============================================================================
// Bethe-Bloch (PDG eq. 34.5), heavy charged particles
// ==============================================================================
//   -dE/dx = K z^2 (Z/A) rho (1/beta^2) [ 1/2 ln(2 m_e c^2 beta^2 gamma^2 Tmax / I^2)
//                                         - beta^2 - delta(beta*gamma)/2 ]
// plus, for spin-1/2 particles (mu, tau, p), the term + (1/8) (Tmax / E)^2 inside
// the bracket (G4MuBetheBlochModel / G4BetheBlochModel have it). It comes from the
// spin-1/2 cross section (the factor 1 - beta^2 T/Tmax + T^2/(2E^2) of PDG eq. 34.6
// integrated up to Tmax); it matters only when Tmax is close to E, i.e. at high
// energy: +0.5 % for a 100 GeV muon in lead.
//
// Not included: shell correction, Barkas (z^3) and Bloch (z^4) terms. They matter
// at low energy only, which is why the PDG table warns "below 10 MeV not dependable"
// for its own values too.
double Ionization::computeBetheBlochDEDX(const Particle& particle, double kineticEnergy, const Material& material) const {
    const double mass = particle.getMass();
    const double charge = particle.getCharge();
    const double totalEnergy = kineticEnergy + mass;
    const double gamma = totalEnergy / mass;
    const double betaGamma2 = kineticEnergy * (kineticEnergy + 2.0 * mass) / (mass * mass);
    const double beta2 = betaGamma2 / (gamma * gamma);

    const double tmax = computeMaxSecondaryEnergy(particle, kineticEnergy);
    const double meanExcitation = material.getMeanExcEn();
    const double delta = material.DensityCorrection(0.5 * std::log10(betaGamma2)); // x = log10(beta*gamma)

    double bracket = 0.5 * std::log(2.0 * m_e * betaGamma2 * tmax / (meanExcitation * meanExcitation))
                   - beta2 - 0.5 * delta;
    if (particle.getSpin() == 0.5) {
        const double ratio = tmax / totalEnergy;
        bracket += 0.125 * ratio * ratio;
    }

    return K * charge * charge * material.getZoverA() * material.getDensity() / beta2 * bracket;
}

// ==============================================================================
// Electrons / positrons: Berger-Seltzer formula (ICRU Report 37)
// ==============================================================================
// Bethe-Bloch does not apply to electrons: the projectile has the same mass as
// the target (so Tmax is T/2 or T, not 2 m_e beta^2 gamma^2), and for e- the two
// are identical particles (Moller scattering, with exchange). The collision
// stopping power then is, with tau = T / m_e:
//
//   -dE/dx = (K/2) (Z/A) rho (1/beta^2) [ ln( tau^2 (tau+2) / (2 (I/m_e)^2) ) + F(tau) - delta ]
//
//   F-(tau) = 1 - beta^2 + [ tau^2/8 - (2 tau + 1) ln2 ] / (tau+1)^2
//   F+(tau) = 2 ln2 - (beta^2/12) [ 23 + 14/(tau+2) + 10/(tau+2)^2 + 4/(tau+2)^3 ]
double Ionization::computeElectronDEDX(const Particle& particle, double kineticEnergy, const Material& material) const {
    const double tau = kineticEnergy / m_e;
    const double gamma = tau + 1.0;
    const double betaGamma2 = tau * (tau + 2.0);
    const double beta2 = betaGamma2 / (gamma * gamma);
    const double ln2 = std::log(2.0);

    double f;
    if (particle.getPdgId() == 11) {
        f = 1.0 - beta2 + (tau * tau / 8.0 - (2.0 * tau + 1.0) * ln2) / (gamma * gamma);
    } else {
        const double y = 1.0 / (tau + 2.0);
        f = 2.0 * ln2 - beta2 / 12.0 * (23.0 + 14.0 * y + 10.0 * y * y + 4.0 * y * y * y);
    }

    const double iOverM = material.getMeanExcEn() / m_e;
    const double delta = material.DensityCorrection(0.5 * std::log10(betaGamma2));
    const double bracket = std::log(tau * tau * (tau + 2.0) / (2.0 * iOverM * iOverM)) + f - delta;

    return 0.5 * K * material.getZoverA() * material.getDensity() / beta2 * bracket;
}

// ==============================================================================
// Range table:  R(T) = integral_0^T dT' / (dE/dx)(T')
// ==============================================================================
// Why a table? The particle needs R(T) at every step (the step limit) and the
// integral is too expensive to redo each time. Geant4 builds the same tables at
// initialisation (that is the "### Run 0 starts" pause: physics tables).
//
// The integral in ln(T):  dT = T d(lnT), so  R = integral (T / S(T)) d(lnT).
// On an equally spaced ln(T) grid the integrand is smooth, and the trapezoid rule
// with 50 points per decade is accurate to ~1e-4.
//
// The piece below the first grid energy Tmin: there S ~ sqrt(T) (see computeDEDX),
// and integral_0^Tmin dT / (S(Tmin) sqrt(T/Tmin)) = 2 Tmin / S(Tmin).
const Ionization::RangeTable& Ionization::getRangeTable(const Particle& particle, const Material& material) {
    const auto key = std::make_pair(&particle, &material);
    auto found = rangeTables.find(key);
    if (found != rangeTables.end()) {
        return found->second;
    }

    RangeTable table;
    const int nBins = static_cast<int>(std::lround(binsPerDecade * std::log10(tableEnergyMax / tableEnergyMin)));
    table.lnEnergyMin = std::log(tableEnergyMin);
    table.dLnEnergy = std::log(tableEnergyMax / tableEnergyMin) / nBins;
    table.range.resize(nBins + 1);
    table.slope.resize(nBins + 1);


    auto slope = [&](double lnEnergy) {
        const double energy = std::exp(lnEnergy);
        return energy / computeContinuousDEDX(particle, energy, material);
    };

    table.slope[0] = slope(table.lnEnergyMin);
    table.range[0] = 2.0 * table.slope[0];

    for (int i = 1; i <= nBins; ++i) {
        const double lnEnergy = table.lnEnergyMin + i * table.dLnEnergy;
        table.slope[i] = slope(lnEnergy);
        const double middle = slope(lnEnergy - 0.5 * table.dLnEnergy);
        table.range[i] = table.range[i - 1] + table.dLnEnergy / 6.0 * (table.slope[i - 1] + 4.0 * middle + table.slope[i]);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Built range table for " << particle.getName()
                  << " in " << material.getName() << " (" << nBins + 1 << " points)\n";
    }
    return rangeTables.emplace(key, std::move(table)).first->second;
}
namespace {

double hermite(double r0, double m0, double r1, double m1, double f) {
    const double f2 = f * f, f3 = f2 * f;
    return (2.0 * f3 - 3.0 * f2 + 1.0) * r0 + (f3 - 2.0 * f2 + f) * m0 + (-2.0 * f3 + 3.0 * f2) * r1 + (f3 - f2) * m1;
}

double hermiteDerivative(double r0, double m0, double r1, double m1, double f) {
    const double f2 = f * f;
    return (6.0 * f2 - 6.0 * f) * r0 + (3.0 * f2 - 4.0 * f + 1.0) * m0 + (-6.0 * f2 + 6.0 * f) * r1 + (3.0 * f2 - 2.0 * f) * m1;
}

} // namespace


// Using Hermite
double Ionization::getRange(const Particle& particle, double kineticEnergy, const Material& material) {
    if (computeContinuousDEDX(particle, tableEnergyMin, material) <= 0.0) {
        return GeoDefs::kInfinity; // vacuum
    }
    const RangeTable& table = getRangeTable(particle, material);

    if (kineticEnergy <= tableEnergyMin) {
        return table.range.front() * std::sqrt(kineticEnergy / tableEnergyMin); // R ~ sqrt(T) there
    }
    if (kineticEnergy >= tableEnergyMax) {
        return table.range.back() + (kineticEnergy - tableEnergyMax) / computeContinuousDEDX(particle, tableEnergyMax, material);
    }

    const double position = (std::log(kineticEnergy) - table.lnEnergyMin) / table.dLnEnergy;
    const int bin = std::min(static_cast<int>(position), static_cast<int>(table.range.size()) - 2);
    const double fraction = position - bin;
    // Hermite, was table.range[bin] + fraction * (table.range[bin + 1] - table.range[bin])
    const double h = table.dLnEnergy;
    return hermite(table.range[bin], h * table.slope[bin], table.range[bin + 1], h * table.slope[bin + 1], fraction);
}

// T(R): the exact inverse of the interpolation above
double Ionization::getKineticEnergyFromRange(const Particle& particle, double range, const Material& material) {
    const RangeTable& table = getRangeTable(particle, material);

    if (range <= table.range.front()) {
        const double ratio = range / table.range.front();
        return tableEnergyMin * ratio * ratio;
    }
    if (range >= table.range.back()) {
        return tableEnergyMax + (range - table.range.back()) * computeContinuousDEDX(particle, tableEnergyMax, material);
    }

    // R grows with T, so the table is sorted: binary search for the bin
    const auto upper = std::upper_bound(table.range.begin(), table.range.end(), range);
    const int bin = static_cast<int>(upper - table.range.begin()) - 1;
    const double h = table.dLnEnergy;
    const double r0 = table.range[bin], r1 = table.range[bin + 1];
    const double m0 = h * table.slope[bin], m1 = h * table.slope[bin + 1];
    double fraction = (range - r0) / (r1 - r0);
    for (int iteration = 0; iteration < 20; ++iteration) {
        const double step = (hermite(r0, m0, r1, m1, fraction) - range) / hermiteDerivative(r0, m0, r1, m1, fraction);
        fraction = std::min(1.0, std::max(0.0, fraction - step));
        if (std::abs(step) < 1e-15) break;
    }
    return std::exp(table.lnEnergyMin + (bin + fraction) * table.dLnEnergy);
}

// ==============================================================================
// Along-step GPIL: how far may the particle go before dE/dx changes too much?
// ==============================================================================
double Ionization::alongStepGetPhysicalInteractionLength(const Track& track, const Material& material) {
    const double range = getRange(*track.getParticleDefinition(), track.getKineticEnergy(), material);
    if (range >= GeoDefs::kInfinity) {
        return GeoDefs::kInfinity; // vacuum: no energy loss, no reason to limit the step
    }

    // Close to the end: let the particle stop in one step (otherwise it would
    // keep losing 20 % of an ever smaller range and take forever).
    if (range <= finalRange) {
        return range;
    }
    return dRmaxRange * range + finalRange * (1.0 - dRmaxRange) * (2.0 - finalRange / range);
}

// ==============================================================================
// Along-step DoIt: continuous energy loss over the step actually taken
// ==============================================================================
void Ionization::alongStepDoIt(Track& track, const Material& material, double stepLength) {
    const Particle& particle = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();
    const double range = getRange(particle, kineticEnergy, material);
    if (range >= GeoDefs::kInfinity) {
        return; // vacuum
    }

    // ------------------------------------------------------------------
    // 1. Mean energy loss of this step
    // ------------------------------------------------------------------
    double energyLoss;
    if (stepLength >= range) {
        energyLoss = kineticEnergy; // the step reaches the end of the range: stop here
    } else {
        // Short step: dE/dx is nearly constant over it
        energyLoss = computeContinuousDEDX(particle, kineticEnergy, material) * stepLength;

        // Long step: dE/dx grows as the particle slows down, so dE/dx * step
        // underestimates the loss. Use the range table instead: after the step the
        // remaining range is R - step, and its energy is T(R - step).
        if (energyLoss > linLossLimit * kineticEnergy) {
            energyLoss = kineticEnergy - getKineticEnergyFromRange(particle, range - stepLength, material);
        }

        // ------------------------------------------------------------------
        // 2. Fluctuation: the actual loss of THIS step around the mean
        // ------------------------------------------------------------------
        // Not for the last step: the particle stops, it loses exactly what it has.
        if (fluctuationModel != FluctuationModel::none && energyLoss + lowestKineticEnergy < kineticEnergy) {
            energyLoss = sampleFluctuations(particle, kineticEnergy, material, energyLoss, stepLength);
        }
    }

    // Tracking cut: if what would be left is tiny, deposit everything here.
    // Without it, rounding can leave e.g. 1e-12 GeV and a few useless steps.
    if (kineticEnergy - energyLoss < lowestKineticEnergy) {
        energyLoss = kineticEnergy;
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Active! "
                  << "Loss: " << energyLoss << " GeV (dE/dx: "
                  << computeContinuousDEDX(particle, kineticEnergy, material) << " GeV/cm)\n";
    }

    track.deductEnergy(energyLoss);
}

// ==============================================================================
// Energy-loss fluctuations
// ==============================================================================
// The mean loss comes from very many soft collisions plus a few rare hard ones
// (up to Tmax). Two numbers decide the shape (PDG section 34.2.9):
//
//   xi    = (K/2) z^2 (Z/A) rho x / beta^2     typical energy of the "few hardest"
//                                              collisions in thickness x
//   kappa = xi / Tmax
//
//   kappa >> 1 (thick layer, slow particle): many collisions of all sizes ->
//                central limit theorem -> GAUSSIAN with Bohr's variance
//                  sigma^2 = xi Tmax (1 - beta^2/2)
//   kappa << 1 (thin layer, fast particle): the loss is dominated by a handful of
//                hard collisions with a 1/T^2 spectrum -> LANDAU: a sharp peak
//                (most probable loss) and a long tail to high losses.
//
// Both samplers below are built so that the MEAN of the sampled loss is exactly
// meanLoss. That is essential: the mean range must not depend on the model.
double Ionization::sampleFluctuations(const Particle& particle, double kineticEnergy, const Material& material,
                                      double meanLoss, double stepLength) const {
    if (fluctuationModel == FluctuationModel::none || meanLoss <= 0.0) {
        return meanLoss;
    }

    const double mass = particle.getMass();
    const double charge = particle.getCharge();
    const double gamma = 1.0 + kineticEnergy / mass;
    const double beta2 = 1.0 - 1.0 / (gamma * gamma);
    const double tmax = std::min(computeMaxSecondaryEnergy(particle, kineticEnergy), kineticEnergy);
    const double tUp = std::min(tmax, deltaRayCut);

    const double xi = 0.5 * K * charge * charge * material.getZoverA() * material.getDensity() * stepLength / beta2;

    // ------------------------------------------------------------------
    // LANDAU
    // ------------------------------------------------------------------
    // Landau's variable: loss = Delta_p + xi (lambda - lambda_p), with lambda_p =
    // -0.22278 the peak of the standard Landau density. The true spectrum ends at
    // Tmax, so lambda is cut at lambda_max = lambda_p + (Tmax - Delta_p) / xi.
    //
    // [INFO] A Landau distribution has NO mean (its tail ~ 1/lambda^2 makes the
    //        integral diverge). Sampled without the cut, the average loss per
    //        step would drift with the number of samples, and ranges would be
    //        nonsense. With the cut, the mean exists, but it is not automatically
    //        equal to the Bethe-Bloch mean. So we shift the distribution:
    //            loss = meanLoss + xi (lambda - <lambda>_cut)
    //        Mean = meanLoss by construction; the shape (width xi, the tail) is
    //        Landau's. The PEAK position is then a prediction, and the validation
    //        program checks it against PDG's Delta_p formula.
    //
    // Landau theory needs lambda_max large (thin layer, kappa << 1). If the layer
    // is too thick for it, use the Gaussian below. (In between, the Vavilov
    // distribution is the right one; Geant4's G4UniversalFluctuation handles all
    // regimes with one more complex model.)
    if (fluctuationModel == FluctuationModel::landau) {
        const double lambdaPeak = -0.22278;
        const double mostProbableLoss = computeMostProbableLoss(particle, kineticEnergy, material, stepLength);
        const double lambdaMax = lambdaPeak + (tUp - mostProbableLoss) / xi;
        if (lambdaMax >= 10.0) {
            const double loss = sampleLandauLoss(meanLoss, xi, lambdaMax);
            // The left side of a Landau is short but not zero: very rarely the
            // shift makes the loss slightly negative. Losing nothing is the limit.
            return std::min(std::max(loss, 0.0), kineticEnergy);
        }
    }

    // ------------------------------------------------------------------
    // GAUSS (Bohr)
    // ------------------------------------------------------------------
    const double sigma2 = xi * tUp * (1.0 - 0.5 * beta2 * tUp / tmax);
    double loss;
    if (meanLoss * meanLoss > 4.0 * sigma2) {
        // mean > 2 sigma: a Gaussian, cut SYMMETRICALLY to [0, 2 mean] so the
        // mean stays meanLoss (cutting only the negative side would raise it).
        const double sigma = std::sqrt(sigma2);
        do {
            loss = RandomEngine::randGauss(meanLoss, sigma);
        } while (loss < 0.0 || loss > 2.0 * meanLoss);
    } else {
        // Wider than that: a cut Gaussian would no longer have the right width.
        // A Gamma distribution with the same mean and variance is positive by
        // construction (G4UniversalFluctuation does the same):
        //   shape n = mean^2 / sigma^2,  loss = mean * Gamma(n, 1) / n
        const double shape = meanLoss * meanLoss / sigma2;
        loss = meanLoss * RandomEngine::randGamma(shape, 1.0) / shape;
    }
    return std::min(loss, kineticEnergy);
}

// ==============================================================================
// Landau most probable loss (PDG eq. 34.12)
// ==============================================================================
//   Delta_p = xi [ ln(2 m_e c^2 beta^2 gamma^2 / I) + ln(xi / I) + j - beta^2 - delta ],  j = 0.200
double Ionization::computeMostProbableLoss(const Particle& particle, double kineticEnergy, const Material& material,
                                           double stepLength) const {
    const double mass = particle.getMass();
    const double charge = particle.getCharge();
    const double betaGamma2 = kineticEnergy * (kineticEnergy + 2.0 * mass) / (mass * mass);
    const double beta2 = betaGamma2 / (1.0 + betaGamma2);
    const double meanExcitation = material.getMeanExcEn();
    const double delta = material.DensityCorrection(0.5 * std::log10(betaGamma2));

    const double xi = 0.5 * K * charge * charge * material.getZoverA() * material.getDensity() * stepLength / beta2;
    return xi * (std::log(2.0 * m_e * betaGamma2 / meanExcitation) + std::log(xi / meanExcitation)
                 + 0.200 - beta2 - delta);
}

// ==============================================================================
// Truncated Landau sampling
// ==============================================================================
namespace {

// Mean of the standard Landau density cut at lambdaMax:
//     <lambda>_cut = integral_{-inf}^{lambdaMax} lambda f(lambda) dlambda / F(lambdaMax)
// It depends only on lambdaMax, so it is tabulated ONCE (first call), on a
// logarithmic grid from 1 to 1e12, and interpolated afterwards.
class TruncatedLandauMean {
private:
    static constexpr int pointsPerDecade = 100;
    static constexpr int decades = 12;
    std::vector<double> mean; // <lambda>_cut at lambdaMax = 10^(i / pointsPerDecade)

    static double integrand(double lambda) { return lambda * ROOT::Math::landau_pdf(lambda); }

public:
    TruncatedLandauMean() {
        // integral of lambda f from -10 (f is ~1e-1000 there: exp(-exp(10))) to 1, Simpson
        const int nLinear = 22000;
        const double lo = -10.0, hi = 1.0, h = (hi - lo) / nLinear;
        double integral = integrand(lo) + integrand(hi);
        for (int i = 1; i < nLinear; ++i) {
            integral += (i % 2 == 1 ? 4.0 : 2.0) * integrand(lo + i * h);
        }
        integral *= h / 3.0;

        // then decade by decade on the log grid; in u = ln(lambda), dlambda = lambda du
        const int nPoints = decades * pointsPerDecade + 1;
        mean.resize(nPoints);
        mean[0] = integral / ROOT::Math::landau_cdf(1.0);
        const double du = std::log(10.0) / pointsPerDecade;
        for (int i = 1; i < nPoints; ++i) {
            const double u0 = (i - 1) * du;
            const int nSub = 8; // Simpson inside each grid cell
            const double hu = du / nSub;
            double cell = 0.0;
            for (int j = 0; j <= nSub; ++j) {
                const double lambda = std::exp(u0 + j * hu);
                const double weight = (j == 0 || j == nSub) ? 1.0 : (j % 2 == 1 ? 4.0 : 2.0);
                cell += weight * integrand(lambda) * lambda;
            }
            integral += cell * hu / 3.0;
            mean[i] = integral / ROOT::Math::landau_cdf(std::exp(i * du));
        }
    }

    double operator()(double lambdaMax) const {
        const double position = std::log10(lambdaMax) * pointsPerDecade;
        if (position <= 0.0) return mean.front();
        if (position >= mean.size() - 1) return mean.back(); // lambdaMax > 1e12: never in practice
        const int bin = static_cast<int>(position);
        const double fraction = position - bin;
        return mean[bin] + fraction * (mean[bin + 1] - mean[bin]);
    }
};

} // namespace

// Draw lambda from the standard Landau cut at lambdaMax (inverse-CDF method:
// u uniform in (0, F(lambdaMax)], lambda = F^-1(u)), then shift it to the mean.
// The uniform number comes from OUR RandomEngine, so the run stays reproducible
// from the one seed (ROOT's gRandom->Landau() would use ROOT's own generator).
double Ionization::sampleLandauLoss(double meanLoss, double xi, double lambdaMax) const {
    static const TruncatedLandauMean truncatedMean; // built at the first call only

    const double cdfMax = ROOT::Math::landau_cdf(lambdaMax);
    const double u = (1.0 - RandomEngine::uniNumGen()) * cdfMax; // in (0, cdfMax]
    const double lambda = ROOT::Math::landau_quantile(u);
    return meanLoss + xi * (lambda - truncatedMean(lambdaMax));
}
// ==============================================================================
// Delta rays
// ==============================================================================
// Notation: T, E = kinetic and total energy of the incident particle, T' = kinetic
// energy given to the electron, eps = T'/T (for e-, e+). Every cross section is
// written as
//     dsigma/dT' = P z(T') / T'^2,     P = (K/2) z^2 (Z/A) rho / beta^2
// P / T'^2 is Rutherford scattering on a free electron ((K/2)(Z/A) rho is
// 2 pi r_e^2 m_e c^2 times the electron density), z(T') ~ 1 is the correction for
// spin, recoil and (for e-) identical particles. Then, in closed form below:
//     sigma(T_cut)  = P integral z / T'^2 dT'     number of delta rays per cm
//     dE/dx_hard    = P integral z / T'  dT'     energy they carry per cm
// integrated from T_cut to Tmax.
namespace {

// Bhabha's coefficients, y = 1 / (gamma + 1) (G4MollerBhabhaModel)
struct BhabhaCoefficients {
    double b1, b2, b3, b4;
    explicit BhabhaCoefficients(double gamma) {
        const double y = 1.0 / (1.0 + gamma);
        const double y12 = 1.0 - 2.0 * y;
        b1 = 2.0 - y * y;
        b2 = y12 * (3.0 + y * y);
        b4 = y12 * y12 * y12;
        b3 = b4 + y12 * y12;
    }
};

} // namespace

double Ionization::computeDeltaRayPrefactor(const Particle& particle, double kineticEnergy, const Material& material) {
    const double mass = particle.getMass();
    const double charge = particle.getCharge();
    const double betaGamma2 = kineticEnergy * (kineticEnergy + 2.0 * mass) / (mass * mass);
    const double beta2 = betaGamma2 / (1.0 + betaGamma2);
    return 0.5 * K * charge * charge * material.getZoverA() * material.getDensity() / beta2;
}

// ------------------------------------------------------------------------------
// The shape z(T')
// ------------------------------------------------------------------------------
//   heavy, spin 0 or 1/2 (PDG eq. 34.6):  z = 1 - beta^2 T'/Tmax  [+ T'^2 / (2 E^2) for spin 1/2]
//   Moller (e- e-), g = (2 gamma - 1) / gamma^2:
//       dsigma/deps = P/T [ (1 - g) + 1/eps^2 - g/eps + 1/(1-eps)^2 - g/(1-eps) ]
//       The 1/(1-eps)^2 terms are exchange: the electrons are identical, so
//       "the incident one kept T'" also counts. That is why eps <= 1/2 only.
//       z = eps^2 [ ... ] = 1 - g eps + eps^2 (1 - g + (1 - g (1-eps)) / (1-eps)^2)
//   Bhabha (e+ e-):
//       dsigma/deps = P beta^2/T [ 1/(beta^2 eps^2) - b1/eps + b2 - b3 eps + b4 eps^2 ]
//       (annihilation diagram included), z = 1 - beta^2 eps (b1 - eps (b2 - eps (b3 - eps b4)))
double Ionization::computeDeltaRayShape(const Particle& particle, double kineticEnergy, double deltaEnergy) const {
    const int pdg = particle.getPdgId();
    if (pdg == 11 || pdg == -11) {
        const double gamma = 1.0 + kineticEnergy / m_e;
        const double eps = deltaEnergy / kineticEnergy;
        if (pdg == 11) {
            const double g = (2.0 * gamma - 1.0) / (gamma * gamma);
            const double y = 1.0 - eps;
            return 1.0 - g * eps + eps * eps * (1.0 - g + (1.0 - g * y) / (y * y));
        }
        const double beta2 = 1.0 - 1.0 / (gamma * gamma);
        const BhabhaCoefficients b(gamma);
        return 1.0 - beta2 * eps * (b.b1 - eps * (b.b2 - eps * (b.b3 - eps * b.b4)));
    }
    const double mass = particle.getMass();
    const double totalEnergy = kineticEnergy + mass;
    const double beta2 = kineticEnergy * (kineticEnergy + 2.0 * mass) / (totalEnergy * totalEnergy);
    double z = 1.0 - beta2 * deltaEnergy / computeMaxSecondaryEnergy(particle, kineticEnergy);
    if (particle.getSpin() == 0.5) {
        const double r = deltaEnergy / totalEnergy;
        z += 0.5 * r * r;
    }
    return z;
}

double Ionization::computeDeltaRayDifferentialCrossSection(const Particle& particle, double kineticEnergy, double deltaEnergy,
                                                           const Material& material) const {
    if (deltaEnergy <= 0.0 || deltaEnergy > getMaxDeltaRayEnergy(particle, kineticEnergy)) {
        return 0.0;
    }
    return computeDeltaRayPrefactor(particle, kineticEnergy, material) * computeDeltaRayShape(particle, kineticEnergy, deltaEnergy)
           / (deltaEnergy * deltaEnergy);
}

// ------------------------------------------------------------------------------
// sigma(T_cut) = P integral z / T'^2 dT'     (x1 = T_cut/T, x2 = Tmax/T)
// ------------------------------------------------------------------------------
//   heavy:  P [ 1/T_cut - 1/Tmax - beta^2 ln(Tmax/T_cut) / Tmax  (+ (Tmax - T_cut) / (2 E^2)) ]
//   Moller: P/T [ (x2 - x1) (1 - g + 1/(x1 x2) + 1/((1-x1)(1-x2))) - g ln(x2 (1-x1) / (x1 (1-x2))) ]
//   Bhabha: P beta^2/T [ (x2 - x1) (1/(beta^2 x1 x2) + b2 - b3 (x1+x2)/2 + b4 (x1^2 + x1 x2 + x2^2)/3)
//                        - b1 ln(x2/x1) ]
// (G4BetheBlochModel / G4MollerBhabhaModel::ComputeCrossSectionPerElectron)
double Ionization::computeDeltaRayCrossSection(const Particle& particle, double kineticEnergy, const Material& material,
                                               double cut) const {
    const double tmax = getMaxDeltaRayEnergy(particle, kineticEnergy);
    if (cut >= tmax) {
        return 0.0;
    }
    const double p = computeDeltaRayPrefactor(particle, kineticEnergy, material);
    const int pdg = particle.getPdgId();
    if (pdg == 11 || pdg == -11) {
        const double gamma = 1.0 + kineticEnergy / m_e;
        const double x1 = cut / kineticEnergy;
        const double x2 = tmax / kineticEnergy;
        if (pdg == 11) {
            const double g = (2.0 * gamma - 1.0) / (gamma * gamma);
            return p / kineticEnergy * ((x2 - x1) * (1.0 - g + 1.0 / (x1 * x2) + 1.0 / ((1.0 - x1) * (1.0 - x2)))
                                        - g * std::log(x2 * (1.0 - x1) / (x1 * (1.0 - x2))));
        }
        const double beta2 = 1.0 - 1.0 / (gamma * gamma);
        const BhabhaCoefficients b(gamma);
        return p * beta2 / kineticEnergy
               * ((x2 - x1) * (1.0 / (beta2 * x1 * x2) + b.b2 - 0.5 * b.b3 * (x1 + x2) + b.b4 * (x1 * x1 + x1 * x2 + x2 * x2) / 3.0)
                  - b.b1 * std::log(x2 / x1));
    }
    const double totalEnergy = kineticEnergy + particle.getMass();
    const double beta2 = kineticEnergy * (kineticEnergy + 2.0 * particle.getMass()) / (totalEnergy * totalEnergy);
    double s = (tmax - cut) / (cut * tmax) - beta2 * std::log(tmax / cut) / computeMaxSecondaryEnergy(particle, kineticEnergy);
    if (particle.getSpin() == 0.5) {
        s += 0.5 * (tmax - cut) / (totalEnergy * totalEnergy);
    }
    return p * s;
}

// ------------------------------------------------------------------------------
// dE/dx_hard = P integral z / T' dT'
// ------------------------------------------------------------------------------
//   heavy:  P [ ln(Tmax/T_cut) - beta^2 (Tmax - T_cut) / Tmax  (+ (Tmax^2 - T_cut^2) / (4 E^2)) ]
//   Moller: P [ F(x2) - F(x1) ],  F(x) = (1 - g) x^2/2 + ln x + 1/(1-x) + (1 + g) ln(1-x)
//   Bhabha: P [ ln(x2/x1) - beta^2 ( b1 (x2-x1) - b2 (x2^2-x1^2)/2 + b3 (x2^3-x1^3)/3 - b4 (x2^4-x1^4)/4 ) ]
// For the heavy particles, unrestricted Bethe-Bloch minus this is exactly the
// restricted Bethe-Bloch of PDG eq. 34.7 (G4BetheBlochModel::ComputeDEDXPerVolume).
// For e- it gives the restricted Berger-Seltzer formula (G4MollerBhabhaModel), the
// validation checks both.
double Ionization::computeDeltaRayEnergyLoss(const Particle& particle, double kineticEnergy, const Material& material,
                                             double cut) const {
    const double tmax = getMaxDeltaRayEnergy(particle, kineticEnergy);
    if (cut >= tmax) {
        return 0.0;
    }
    const double p = computeDeltaRayPrefactor(particle, kineticEnergy, material);
    const int pdg = particle.getPdgId();
    if (pdg == 11 || pdg == -11) {
        const double gamma = 1.0 + kineticEnergy / m_e;
        const double x1 = cut / kineticEnergy;
        const double x2 = tmax / kineticEnergy;
        if (pdg == 11) {
            const double g = (2.0 * gamma - 1.0) / (gamma * gamma);
            auto f = [g](double x) { return 0.5 * (1.0 - g) * x * x + std::log(x) + 1.0 / (1.0 - x) + (1.0 + g) * std::log(1.0 - x); };
            return p * (f(x2) - f(x1));
        }
        const double beta2 = 1.0 - 1.0 / (gamma * gamma);
        const BhabhaCoefficients b(gamma);
        auto power = [x1, x2](int n) { return std::pow(x2, n) - std::pow(x1, n); };
        return p * (std::log(x2 / x1)
                    - beta2 * (b.b1 * power(1) - 0.5 * b.b2 * power(2) + b.b3 * power(3) / 3.0 - 0.25 * b.b4 * power(4)));
    }
    const double totalEnergy = kineticEnergy + particle.getMass();
    const double beta2 = kineticEnergy * (kineticEnergy + 2.0 * particle.getMass()) / (totalEnergy * totalEnergy);
    double s = std::log(tmax / cut) - beta2 * (tmax - cut) / computeMaxSecondaryEnergy(particle, kineticEnergy);
    if (particle.getSpin() == 0.5) {
        s += 0.25 * (tmax * tmax - cut * cut) / (totalEnergy * totalEnergy);
    }
    return p * s;
}

// Restricted = unrestricted - what the delta rays take. Written like this (and not
// as a separate formula) so that the two halves of the process add up EXACTLY to
// the stopping power that was validated against PDG and ESTAR in stage 2.
double Ionization::computeRestrictedDEDX(const Particle& particle, double kineticEnergy, const Material& material,
                                         double cut) const {
    const double full = computeDEDX(particle, kineticEnergy, material);
    if (full <= 0.0) {
        return 0.0;
    }
    return std::max(0.0, full - computeDeltaRayEnergyLoss(particle, kineticEnergy, material, cut));
}

// ------------------------------------------------------------------------------
// Sampling T'
// ------------------------------------------------------------------------------
// The 1/T'^2 part by inversion (1/T' uniform between 1/Tmax and 1/T_cut):
//     T' = T_cut Tmax / (T_cut (1 - u) + Tmax u)
// then accept with z(T') / zMax, where zMax is an upper bound of z:
//   heavy:  1 (+ Tmax^2 / (2 E^2) for spin 1/2)
//   Moller: the larger of z at the two ends. z has no maximum inside the interval
//           (checked numerically for gamma = 1.0001 ... 1e5, every cut)
//   Bhabha: each power of eps taken at the end that makes it largest
// (as G4BetheBlochModel / G4MollerBhabhaModel::SampleSecondaries do)
double Ionization::sampleDeltaRayEnergy(const Particle& particle, double kineticEnergy, double cut) const {
    const double tmax = getMaxDeltaRayEnergy(particle, kineticEnergy);
    if (cut >= tmax) {
        return 0.0;
    }
    const int pdg = particle.getPdgId();
    double zMax;
    if (pdg == 11) {
        zMax = std::max(computeDeltaRayShape(particle, kineticEnergy, cut), computeDeltaRayShape(particle, kineticEnergy, tmax));
    } else if (pdg == -11) {
        const double gamma = 1.0 + kineticEnergy / m_e;
        const double beta2 = 1.0 - 1.0 / (gamma * gamma);
        const BhabhaCoefficients b(gamma);
        const double x1 = cut / kineticEnergy;
        const double x2 = tmax / kineticEnergy;
        zMax = 1.0 + beta2 * (b.b4 * x2 * x2 * x2 * x2 - b.b3 * x1 * x1 * x1 + b.b2 * x2 * x2 - b.b1 * x1);
    } else {
        zMax = 1.0;
        if (particle.getSpin() == 0.5) {
            const double r = tmax / (kineticEnergy + particle.getMass());
            zMax += 0.5 * r * r;
        }
    }

    double deltaEnergy;
    do {
        const double u = RandomEngine::uniNumGen();
        deltaEnergy = cut * tmax / (cut * (1.0 - u) + tmax * u);
    } while (zMax * RandomEngine::uniNumGen() > computeDeltaRayShape(particle, kineticEnergy, deltaEnergy));
    return deltaEnergy;
}

// ------------------------------------------------------------------------------
// Post-step GPIL: the integral approach
// ------------------------------------------------------------------------------
// The particle loses energy DURING the step, so sigma(T) changes along it. Taking
// sigma at the pre-step energy (what every other discrete process here does) is
// biased, and the bias has a sign:
//   cut 1 MeV, mu 200 MeV: sigma FALLS as the muon slows (it is 0 below ~42 MeV,
//       where Tmax < cut): pre-step sigma -> too many delta rays -> range 0.3 % short
//   cut 0.1 MeV: sigma RISES as the muon slows (1/beta^2): too few -> range too long
// (found by test 5 of the validation: -4 and +1.5 standard deviations).
//
// Geant4's fix (G4VEnergyLossProcess, "integral approach"), a null-collision method:
//   1. GPIL: sample the distance with sigmaMax >= sigma(T) for every T the step can
//      reach: from the pre-step energy down to the energy after the longest
//      continuous step allowed
//   2. DoIt: at the end of the step, a real delta ray with probability
//      sigma(T_post) / sigmaMax, otherwise nothing happens ("fictitious" interaction)
// This is EXACT for any sigma(T) below sigmaMax: rate of fictitious points times
// acceptance = sigmaMax * sigma(T)/sigmaMax = sigma(T) at every point of the path.
// Scan up in energy in steps of 0.2 % from the first energy with sigma > 0 until
// sigma falls: that is the peak (to 0.2 %, i.e. sigma to ~1e-6 there).
double Ionization::getSigmaMaxEnergies(const Particle& particle, const Material& material) {
    const auto key = std::make_pair(&particle, &material);
    const auto found = sigmaMaxEnergies.find(key);
    if (found != sigmaMaxEnergies.end()) {
        return found->second;
    }
    const double factor = 1.002;
    double energy = tableEnergyMin;
    double previous = computeDeltaRayCrossSection(particle, energy, material, deltaRayCut);
    double peak = GeoDefs::kInfinity;
    while (energy < tableEnergyMax) {
        const double next = computeDeltaRayCrossSection(particle, energy * factor, material, deltaRayCut);
        if (previous > 0.0 && next < previous) {
            peak = energy;
            break;
        }
        previous = next;
        energy *= factor;
    }
    sigmaMaxEnergies.emplace(key, peak);
    return peak;
}

double Ionization::postStepGetPhysicalInteractionLength(const Track& track, const Material& material) {
    if (deltaElectronDef == nullptr) {
        return GeoDefs::kInfinity;
    }
    const Particle& particle = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();

    // Lowest energy the continuous loss of this step can reach (on average): after
    // the longest step the along-step limit allows
    double lowEnergy = 1e-3 * kineticEnergy;
    const double range = getRange(particle, kineticEnergy, material);
    if (range < GeoDefs::kInfinity) {
        const double longestStep = alongStepGetPhysicalInteractionLength(track, material);
        if (longestStep < range) {
            lowEnergy = std::max(lowEnergy, getKineticEnergyFromRange(particle, range - longestStep, material));
        }
    } else {
        lowEnergy = kineticEnergy; // vacuum: no loss, sigma = 0 anyway
    }

    // The largest sigma in [lowEnergy, T]. sigma(T) is 0 below the threshold
    // (Tmax = cut), rises to ONE maximum, then changes slowly (falls like 1/beta^2,
    // rises again a little at high energy). So the largest value in an interval is
    // at one of its ends, or at the peak if the peak is inside.
    sigmaMax = std::max(computeDeltaRayCrossSection(particle, kineticEnergy, material, deltaRayCut),
                        computeDeltaRayCrossSection(particle, lowEnergy, material, deltaRayCut));
    const double peak = getSigmaMaxEnergies(particle, material);
    if (peak > lowEnergy && peak < kineticEnergy) {
        sigmaMax = std::max(sigmaMax, computeDeltaRayCrossSection(particle, peak, material, deltaRayCut));
    }
    // getMeanFreePath (header) returns 1 / sigmaMax
    return DiscreteProcess::postStepGetPhysicalInteractionLength(track, material);
}

// Kinematics of a collision with a free electron at rest (exact, any projectile):
//     cos(theta') = T' (E + m_e) / (p' p)          (G4BetheBlochModel::SampleSecondaries)
// The projectile then has T - T' and momentum p - p' (vectors). With the free-
// electron kinematics |p - p'| is exactly the momentum of T - T', which the
// validation checks.
void Ionization::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>& secondaries) {
    clearNumberOfInteractionLengthLeft();
    const Particle& particle = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();

    // Integral approach, step 2: real delta ray with probability sigma(T_post) / sigmaMax.
    // (Also covers "the along-step loss brought Tmax below the cut": sigma = 0 there.)
    const double sigma = computeDeltaRayCrossSection(particle, kineticEnergy, material, deltaRayCut);
    if (sigma > sigmaMax) {
        ++sigmaMaxViolations;
    }
    if (sigma <= 0.0 || sigma < sigmaMax * RandomEngine::uniNumGen()) {
        return; // fictitious interaction: nothing happens
    }
    const double deltaEnergy = sampleDeltaRayEnergy(particle, kineticEnergy, deltaRayCut);

    const double momentum = track.getMomentumMag();
    const double deltaMomentum = std::sqrt(deltaEnergy * (deltaEnergy + 2.0 * m_e));
    const double cosTheta = std::min(1.0, deltaEnergy * (track.getTotalEnergy() + m_e) / (deltaMomentum * momentum));
    const double sinTheta = std::sqrt((1.0 - cosTheta) * (1.0 + cosTheta));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    Vector3D deltaDirection(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    deltaDirection.rotateUz(track.getDirection());
    secondaries.emplace_back(deltaElectronDef, track.getPosition(), deltaDirection, deltaEnergy);

    // The projectile: momentum conservation. An e+ can give ALL its energy
    // (eps = 1): it then stops (trackStopButAlive, annihilates at rest) and has no
    // direction any more.
    const Vector3D newMomentum = track.getDirection() * momentum - deltaDirection * deltaMomentum;
    track.deductEnergy(deltaEnergy);
    if (track.getKineticEnergy() > 0.0) {
        track.setDirection(newMomentum);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] delta ray of " << deltaEnergy / Units::MeV << " MeV (cos theta "
                  << cosTheta << "), " << particle.getName() << " " << kineticEnergy / Units::MeV << " -> "
                  << track.getKineticEnergy() / Units::MeV << " MeV\n";
    }
}

