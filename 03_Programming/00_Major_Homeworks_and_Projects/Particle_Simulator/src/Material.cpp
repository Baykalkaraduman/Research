#include "Material.h"

#include <cmath>

// Physics: a fast particle polarises the medium, and the polarisation screens
// its far field. Distant collisions get weaker, so dE/dx rises more slowly than
// the bare ln(beta*gamma) of the Bethe formula. delta/2 is subtracted from the
// bracket. Without it, dE/dx in lead at 100 GeV would come out ~25 % too high.
//
//   x = log10(beta*gamma)
//   x >= x1       : delta = 2 ln10 x - cBar
//   x0 <= x < x1  : delta = 2 ln10 x - cBar + a (x1 - x)^k
//   x < x0        : delta = delta0 * 10^(2 (x - x0))     (conductors; 0 for insulators)
double Material::DensityCorrection(double x) const {
    const DensityEffectData& d = densityEffect;
    const double twoLn10 = 2.0 * std::log(10.0);

    if (x >= d.x1) {
        return twoLn10 * x - d.cBar;
    }
    if (x >= d.x0) {
        return twoLn10 * x - d.cBar + d.a * std::pow(d.x1 - x, d.k);
    }
    return d.delta0 * std::pow(10.0, 2.0 * (x - d.x0));
}

double Material::computeCoulombCorrection() const {
    const double az = PhysicalConstants::alpha * Z;
    const double az2 = az * az;
    return az2 * (1.0 / (1.0 + az2) + 0.20206 - 0.0369 * az2 + 0.0083 * az2 * az2 - 0.002 * az2 * az2 * az2);
}

