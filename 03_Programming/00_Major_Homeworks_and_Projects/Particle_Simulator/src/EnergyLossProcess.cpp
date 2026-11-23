#include "EnergyLossProcess.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// ==============================================================================
// Numerical Integration in ln(eps)
// ==============================================================================
// Every dsigma/deps here behaves roughly like 1/eps (brems, photonuclear) or
// falls faster (pairs), so the integrand is smooth in ln(eps):
//     integral g(eps) deps = integral g(eps) eps d(ln eps)
//
// 8-point Gauss-Legendre on sub-intervals of at most 0.25 in ln(eps) (~9 per
// decade). Discontinuities (e.g., k_max' in bremsstrahlung) are split via
// integrateWithBreakpoints so no step jump falls inside a Gauss sub-interval.
namespace {

const double gaussX[8] = {
    0.0198550717512320, 0.1016667612931865, 0.2372337950418355, 0.4082826787521750,
    0.5917173212478250, 0.7627662049581645, 0.8983332387068135, 0.9801449282487680
};

const double gaussW[8] = {
    0.0506142681451880, 0.1111905172266870, 0.1568533229389435, 0.1813418916891810,
    0.1813418916891810, 0.1568533229389435, 0.1111905172266870, 0.0506142681451880
};

// Computes integral_lo^hi g(eps) deps, with g given as a function of eps.
template <typename G>
double integrateInLog(G g, double lo, double hi) {
    if (hi <= lo || lo <= 0.0) {
        return 0.0;
    }

    const double a = std::log(lo);
    const double b = std::log(hi);
    const int n = std::max(4, static_cast<int>(std::ceil((b - a) / 0.25)));
    const double h = (b - a) / n;

    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 8; ++j) {
            const double eps = std::exp(a + (i + gaussX[j]) * h);
            sum += gaussW[j] * g(eps) * eps;
        }
    }

    return sum * h;
}

// Splits [lo, hi] at any interior breakpoints and integrates each piece cleanly.
template <typename G>
double integrateWithBreakpoints(G g, double lo, double hi, std::vector<double> breakpoints) {
    if (hi <= lo || lo <= 0.0) {
        return 0.0;
    }

    std::sort(breakpoints.begin(), breakpoints.end());

    double sum = 0.0;
    double currentLo = lo;
    for (double bp : breakpoints) {
        if (bp > currentLo && bp < hi) {
            sum += integrateInLog(g, currentLo, bp);
            currentLo = bp;
        }
    }
    sum += integrateInLog(g, currentLo, hi);

    return sum;
}

} // namespace

// ==============================================================================
// Internal Helpers: Tables & Interpolation
// ==============================================================================
// 20 points per decade from the lowest energy of the process to 1 PeV. sigma
// and dE/dx_restr change slowly (like ln T) above threshold, so linear
// interpolation in ln(T) is accurate to ~1e-3 there; just above the thresholds
// (pairs: 0.85 GeV) it is coarser, but the cross sections are tiny there.
const EnergyLossProcess::LossTable& EnergyLossProcess::getTable(const Particle& particle, const Material& material) const {
    const TableKey key = std::make_pair(&particle, &material);
    const auto found = tables.find(key);
    if (found != tables.end()) {
        return found->second;
    }

    LossTable table;
    const double energyMin = getLowestKineticEnergy(particle);
    const int nBins = std::max(1, static_cast<int>(std::lround(binsPerDecade * std::log10(tableEnergyMax / energyMin))));

    table.lnEnergyMin = std::log(energyMin);
    table.dLnEnergy = std::log(tableEnergyMax / energyMin) / nBins;
    table.sigma.resize(nBins + 1);
    table.restrictedDEDX.resize(nBins + 1);
    table.majorant.resize(nBins + 1);

    const double n = material.getAtomDensity();
    for (int i = 0; i <= nBins; ++i) {
        const double energy = std::exp(table.lnEnergyMin + i * table.dLnEnergy);
        const double lo = std::max(productionCut, getMinEnergyTransfer(particle, material));
        const double hi = getMaxEnergyTransfer(particle, energy, material);

        table.sigma[i] = n * computeCrossSectionPerAtom(particle, energy, productionCut, energy, material);
        table.restrictedDEDX[i] = n * computeEnergyLossPerAtom(particle, energy, 0.0, productionCut, material);
        table.majorant[i] = computeMajorant(particle, energy, lo, hi, material);
    }

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] Built tables for " << particle.getName()
                  << " in " << material.getName() << " (cut " << productionCut / Units::MeV
                  << " MeV, " << nBins + 1 << " points)\n";
    }

    return tables.emplace(key, std::move(table)).first->second;
}

double EnergyLossProcess::interpolate(const LossTable& table, const std::vector<double>& values, double energy) {
    const double position = (std::log(energy) - table.lnEnergyMin) / table.dLnEnergy;
    if (position <= 0.0) {
        return values.front();
    }
    if (position >= static_cast<double>(values.size() - 1)) {
        return values.back();
    }

    const int bin = static_cast<int>(position);
    const double fraction = position - bin;
    return values[bin] + fraction * (values[bin + 1] - values[bin]);
}

// Must be set before the run: lookup tables depend on the production cut.
void EnergyLossProcess::setProductionCut(double cut) {
    productionCut = cut;
    tables.clear();
}

// ==============================================================================
// Internal Helpers: Majorant & Rejection Sampling
// ==============================================================================
// Envelope 1/eps: ln(eps) uniform in [ln lo, ln hi]. Accept with f(eps) / fMax,
// where f = eps * dsigma/deps is the density in ln(eps). Exact as long as
// fMax >= f everywhere (counted in majorantViolations).
//
// fMax is found by scanning f at 64 points equally spaced in ln(eps), times 1.2.
// A scan is required because pair production rises from the 4 m_e threshold
// before falling, unlike bremsstrahlung which falls monotonically.
double EnergyLossProcess::computeMajorant(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const {
    if (hi <= lo) {
        return 0.0;
    }

    const int n = 64;
    const double a = std::log(lo);
    const double h = (std::log(hi) - a) / n;

    double fMax = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double eps = (i == n) ? hi : std::exp(a + i * h);
        fMax = std::max(fMax, eps * computeDifferentialCrossSectionPerAtom(particle, kineticEnergy, eps, material));
    }

    return 1.2 * fMax;
}

double EnergyLossProcess::sampleWithMajorant(const Particle& particle, double kineticEnergy, double lo, double hi, double fMax, const Material& material) const {
    const double lnRatio = std::log(hi / lo);
    double eps = 0.0;
    double f = 0.0;

    do {
        eps = lo * std::exp(lnRatio * RandomEngine::uniNumGen());
        f = eps * computeDifferentialCrossSectionPerAtom(particle, kineticEnergy, eps, material);
        if (f > fMax) {
            ++majorantViolations;
        }
    } while (f < fMax * RandomEngine::uniNumGen());

    return eps;
}

// ==============================================================================
// DiscreteProcess Overrides
// ==============================================================================
double EnergyLossProcess::getMeanFreePath(const Track& track, const Material& material) const {
    const Particle& particle = *track.getParticleDefinition();
    if (material.getZ() < 1.0 || track.getKineticEnergy() <= getLowestKineticEnergy(particle)) {
        return GeoDefs::kInfinity;
    }

    const LossTable& table = getTable(particle, material);
    const double sigma = interpolate(table, table.sigma, track.getKineticEnergy());
    return sigma > 0.0 ? 1.0 / sigma : GeoDefs::kInfinity;
}

// ==============================================================================
// Stepping Interface: Along-Step (Continuous Loss Below Cut)
// ==============================================================================
// Runs AFTER Ionization in the stepping loop, so alongStepDoIt sees the track
// energy after ionization loss. Using dE/dx at that post-ionization energy alone
// underestimates the loss (~0.5% for 10 GeV muons in 1 m lead with 1 GeV cuts).
// The trapezoid rule between dE/dx(T_pre) (saved here) and dE/dx(T_now) is exact
// for dE/dx linear in T, which radiative losses nearly are (dE/dx ~ b * E).
double EnergyLossProcess::alongStepGetPhysicalInteractionLength(const Track& track, const Material& material) {
    const Particle& particle = *track.getParticleDefinition();
    preStepDEDX = 0.0;

    if (material.getZ() >= 1.0 && track.getKineticEnergy() > getLowestKineticEnergy(particle)) {
        const LossTable& table = getTable(particle, material);
        preStepDEDX = interpolate(table, table.restrictedDEDX, track.getKineticEnergy());
    }

    return GeoDefs::kInfinity;
}

void EnergyLossProcess::alongStepDoIt(Track& track, const Material& material, double stepLength) {
    const Particle& particle = *track.getParticleDefinition();
    if (material.getZ() < 1.0 || track.getKineticEnergy() <= getLowestKineticEnergy(particle)) {
        return;
    }

    const LossTable& table = getTable(particle, material);
    const double dedxNow = interpolate(table, table.restrictedDEDX, track.getKineticEnergy());
    const double loss = std::min(0.5 * (preStepDEDX + dedxNow) * stepLength, track.getKineticEnergy());

    if (loss > 0.0) {
        track.deductEnergy(loss);
    }
}

// ==============================================================================
// Stepping Interface: Post-Step (Discrete Transfer Above Cut)
// ==============================================================================
// 1. Integral approach: sigma(T) grows monotonically with T, so the pre-step
//    sigma (1 / currentInteractionLength) bounds the step. Accept a real
//    interaction with probability sigma(T_post) / sigma(T_pre); otherwise
//    perform a fictitious interaction.
// 2. Majorant from table: take the max of the two surrounding nodes. Since
//    eps * dsigma/deps and [cut, eps_max(T)] grow with T, the upper node bounds
//    every T in its bin.
void EnergyLossProcess::postStepDoIt(Track& track, const Material& material, Step& step, std::vector<Track>& secondaries) {
    clearNumberOfInteractionLengthLeft();

    const Particle& particle = *track.getParticleDefinition();
    const double kineticEnergy = track.getKineticEnergy();
    if (kineticEnergy <= getLowestKineticEnergy(particle)) {
        return; // along-step loss brought it below the model's limit
    }

    const double lo = std::max(productionCut, getMinEnergyTransfer(particle, material));
    const double hi = getMaxEnergyTransfer(particle, kineticEnergy, material);
    if (lo >= hi) {
        return;
    }

    const LossTable& table = getTable(particle, material);
    const double sigmaPost = interpolate(table, table.sigma, kineticEnergy);
    const double sigmaPre = 1.0 / currentInteractionLength;

    if (sigmaPost > sigmaPre * (1.0 + 1e-12)) {
        ++sigmaMaxViolations; // sigma rose during the step: should never happen
    }
    if (sigmaPost < sigmaPre * RandomEngine::uniNumGen()) {
        return; // fictitious interaction
    }

    const double position = (std::log(kineticEnergy) - table.lnEnergyMin) / table.dLnEnergy;
    const int last = static_cast<int>(table.majorant.size()) - 1;
    const int bin = std::min(std::max(static_cast<int>(position), 0), last - 1);

    double fMax = std::max(table.majorant[bin], table.majorant[bin + 1]);
    if (position > last) {
        fMax = computeMajorant(particle, kineticEnergy, lo, hi, material); // above 1 PeV: no table
    }

    const double transfer = sampleWithMajorant(particle, kineticEnergy, lo, hi, fMax, material);
    produceSecondaries(track, material, step, secondaries, transfer);
}

// ==============================================================================
// Physics & Validation Interface
// ==============================================================================
double EnergyLossProcess::computeCrossSectionPerAtom(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const {
    lo = std::max(lo, getMinEnergyTransfer(particle, material));
    hi = std::min(hi, getMaxEnergyTransfer(particle, kineticEnergy, material));
    const auto breakpoints = getIntegrationBreakpoints(particle, kineticEnergy, material);

    return integrateWithBreakpoints(
        [&](double eps) {return computeDifferentialCrossSectionPerAtom(particle, kineticEnergy, eps, material);},
        lo, hi, breakpoints
    );
}

double EnergyLossProcess::computeEnergyLossPerAtom(const Particle& particle, double kineticEnergy, double lo, double hi, const Material& material) const {
    lo = std::max(lo, getMinEnergyTransfer(particle, material));
    hi = std::min(hi, getMaxEnergyTransfer(particle, kineticEnergy, material));
    const auto breakpoints = getIntegrationBreakpoints(particle, kineticEnergy, material);

    return integrateWithBreakpoints(
        [&](double eps) {return eps * computeDifferentialCrossSectionPerAtom(particle, kineticEnergy, eps, material);},
        lo, hi, breakpoints
    );
}

double EnergyLossProcess::computeDEDX(const Particle& particle, double kineticEnergy, const Material& material) const {
    if (material.getZ() < 1.0 || kineticEnergy <= getLowestKineticEnergy(particle)) {
        return 0.0;
    }

    return material.getAtomDensity() * computeEnergyLossPerAtom(particle, kineticEnergy, 0.0, kineticEnergy, material);
}

double EnergyLossProcess::computeRestrictedDEDX(const Particle& particle, double kineticEnergy, const Material& material, double cut) const {
    if (material.getZ() < 1.0 || kineticEnergy <= getLowestKineticEnergy(particle)) {
        return 0.0;
    }

    return material.getAtomDensity() * computeEnergyLossPerAtom(particle, kineticEnergy, 0.0, cut, material);
}

double EnergyLossProcess::sampleEnergyTransfer(const Particle& particle, double kineticEnergy, double cut, const Material& material) const {
    const double lo = std::max(cut, getMinEnergyTransfer(particle, material));
    const double hi = getMaxEnergyTransfer(particle, kineticEnergy, material);
    if (lo >= hi) {
        return 0.0;
    }

    const double fMax = computeMajorant(particle, kineticEnergy, lo, hi, material);
    return sampleWithMajorant(particle, kineticEnergy, lo, hi, fMax, material);
}