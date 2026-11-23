#include "ModifiedMephi.h"
#include "RandomEngine.h"

#include <algorithm>
#include <cmath>

// SampleCosTheta
double ModifiedMephi::sampleCosTheta(double kineticEnergy, double secondaryEnergy, double mass) {
    const double gamma = 1.0 + kineticEnergy / mass;
    const double rMax = gamma * 0.5 * M_PI * std::min(1.0, gamma * mass / secondaryEnergy - 1.0);
    const double rMax2 = rMax * rMax;
    const double x = RandomEngine::uniNumGen() * rMax2 / (1.0 + rMax2);
    
    return std::cos(std::sqrt(x / (1.0 - x)) / gamma);
}
