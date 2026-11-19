#include "ModifiedTsai.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <cmath>

double ModifiedTsai::sampleCosTheta(double kineticEnergy) {
    const double uMax = 2.0 * (1.0 + kineticEnergy / PhysicalConstants::m_e);
    constexpr double a1 = 1.6;
    constexpr double a2 = a1 / 3.0;
    constexpr double border = 0.25;

    double u;
    do {
        const double uu = -std::log((1.0 - RandomEngine::uniNumGen()) * (1.0 - RandomEngine::uniNumGen()));
        u = (border > RandomEngine::uniNumGen()) ? uu * a1 : uu * a2;
    } while (u > uMax);

    return 1.0 - 2.0 * u * u / (uMax * uMax);
}
