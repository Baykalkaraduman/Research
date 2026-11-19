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

    double previous = tableEnergyMin / computeDEDX(particle, tableEnergyMin, material); // T / S(T)
    table.range[0] = 2.0 * previous;
    for (int i = 1; i <= nBins; ++i) {
        const double energy = std::exp(table.lnEnergyMin + i * table.dLnEnergy);
        const double current = energy / computeDEDX(particle, energy, material);
        table.range[i] = table.range[i - 1] + 0.5 * (previous + current) * table.dLnEnergy;
        previous = current;
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Built range table for " << particle.getName()
                  << " in " << material.getName() << " (" << nBins + 1 << " points)\n";
    }
    return rangeTables.emplace(key, std::move(table)).first->second;
}

// R(T): linear interpolation in ln(T)
double Ionization::getRange(const Particle& particle, double kineticEnergy, const Material& material) {
    if (computeDEDX(particle, tableEnergyMin, material) <= 0.0) {
        return GeoDefs::kInfinity; // vacuum
    }
    const RangeTable& table = getRangeTable(particle, material);

    if (kineticEnergy <= tableEnergyMin) {
        return table.range.front() * std::sqrt(kineticEnergy / tableEnergyMin); // R ~ sqrt(T) there
    }
    if (kineticEnergy >= tableEnergyMax) {
        return table.range.back() + (kineticEnergy - tableEnergyMax) / computeDEDX(particle, tableEnergyMax, material);
    }

    const double position = (std::log(kineticEnergy) - table.lnEnergyMin) / table.dLnEnergy;
    const int bin = std::min(static_cast<int>(position), static_cast<int>(table.range.size()) - 2);
    const double fraction = position - bin;
    return table.range[bin] + fraction * (table.range[bin + 1] - table.range[bin]);
}

// T(R): the exact inverse of the interpolation above
double Ionization::getKineticEnergyFromRange(const Particle& particle, double range, const Material& material) {
    const RangeTable& table = getRangeTable(particle, material);

    if (range <= table.range.front()) {
        const double ratio = range / table.range.front();
        return tableEnergyMin * ratio * ratio;
    }
    if (range >= table.range.back()) {
        return tableEnergyMax + (range - table.range.back()) * computeDEDX(particle, tableEnergyMax, material);
    }

    // R grows with T, so the table is sorted: binary search for the bin
    const auto upper = std::upper_bound(table.range.begin(), table.range.end(), range);
    const int bin = static_cast<int>(upper - table.range.begin()) - 1;
    const double fraction = (range - table.range[bin]) / (table.range[bin + 1] - table.range[bin]);
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
        energyLoss = computeDEDX(particle, kineticEnergy, material) * stepLength;

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
                  << computeDEDX(particle, kineticEnergy, material) << " GeV/cm)\n";
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
        const double lambdaMax = lambdaPeak + (tmax - mostProbableLoss) / xi;
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
    const double sigma2 = xi * tmax * (1.0 - 0.5 * beta2);
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
