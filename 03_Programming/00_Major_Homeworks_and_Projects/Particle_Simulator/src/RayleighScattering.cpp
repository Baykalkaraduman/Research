#include "RayleighScattering.h"
#include "PhotoElectricShellData.h"
#include "RandomEngine.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

// ==============================================================================
// Element CDF Cache (Cumulative F^2 in x^2)
// ==============================================================================
namespace {

constexpr double kHcEvAngstrom = 12398.419843320026; // h*c [eV * Angstrom]

struct ElementRayleighCdf {
    // Energy-independent cumulative table A(u) = integral_0^u F^2(sqrt(u')) du'
    // with u = x^2 [1/Angstrom^2]
    std::vector<double> uGrid;
    std::vector<double> f2Grid;
    std::vector<double> cdfGrid;

    explicit ElementRayleighCdf(int z) {
        const PhotoElectricShellData& shellData = PhotoElectricShellData::get(z);

        constexpr int kPointsPerDecade = 100;
        constexpr int kDecades = 10;
        constexpr int kLogPoints = kPointsPerDecade * kDecades + 1;
        constexpr double kMinLogX = -4.0;

        uGrid.reserve(static_cast<std::size_t>(kLogPoints + 1));
        f2Grid.reserve(static_cast<std::size_t>(kLogPoints + 1));
        cdfGrid.reserve(static_cast<std::size_t>(kLogPoints + 1));

        const double f0 = shellData.getFormFactor(0.0);
        uGrid.push_back(0.0);
        f2Grid.push_back(f0 * f0);
        cdfGrid.push_back(0.0);

        for (int i = 0; i < kLogPoints; ++i) {
            const double x = std::pow(10.0, kMinLogX + static_cast<double>(i) / kPointsPerDecade);
            const double u = x * x;
            const double f = shellData.getFormFactor(x);
            const double f2 = f * f;

            const double du = u - uGrid.back();
            const double area = 0.5 * (f2Grid.back() + f2) * du;

            uGrid.push_back(u);
            f2Grid.push_back(f2);
            cdfGrid.push_back(cdfGrid.back() + area);
        }
    }

    [[nodiscard]] double evaluateCdf(double uQuery) const {
        if (uQuery <= 0.0 || uGrid.size() < 2) {
            return 0.0;
        }
        if (uQuery >= uGrid.back()) {
            return cdfGrid.back();
        }

        const auto upper = std::upper_bound(uGrid.begin(), uGrid.end(), uQuery);
        const std::size_t i = static_cast<std::size_t>(upper - uGrid.begin()) - 1;
        const double du = uGrid[i + 1] - uGrid[i];
        const double t = uQuery - uGrid[i];
        const double f2AtQuery = f2Grid[i] + (f2Grid[i + 1] - f2Grid[i]) * (t / du);
        return cdfGrid[i] + 0.5 * (f2Grid[i] + f2AtQuery) * t;
    }

    [[nodiscard]] double sampleX2(double targetA, double uMax) const {
        if (targetA <= 0.0 || uMax <= 0.0) {
            return 0.0;
        }

        const auto upper = std::upper_bound(cdfGrid.begin(), cdfGrid.end(), targetA);
        if (upper == cdfGrid.end()) {
            return std::min(uGrid.back(), uMax);
        }

        const std::size_t i = static_cast<std::size_t>(upper - cdfGrid.begin()) - 1;
        const double du = uGrid[i + 1] - uGrid[i];
        const double f0 = f2Grid[i];
        const double f1 = f2Grid[i + 1];
        const double deltaA = targetA - cdfGrid[i];
        const double slope = (f1 - f0) / du;
        const double disc = std::max(0.0, f0 * f0 + 2.0 * slope * deltaA);
        const double denom = f0 + std::sqrt(disc);
        const double t = (denom > 0.0) ? (2.0 * deltaA) / denom : 0.0;

        return std::clamp(uGrid[i] + t, 0.0, uMax);
    }

    [[nodiscard]] static const ElementRayleighCdf& get(int z) {
        if (z >= 1 && z <= 100) {
            static std::array<std::unique_ptr<ElementRayleighCdf>, 101> fastCache;
            auto& entry = fastCache[static_cast<std::size_t>(z)];
            if (!entry) {
                entry = std::make_unique<ElementRayleighCdf>(z);
            }
            return *entry;
        }

        static std::map<int, std::unique_ptr<ElementRayleighCdf>> extraCache;
        auto& entry = extraCache[z];
        if (!entry) {
            entry = std::make_unique<ElementRayleighCdf>(z);
        }
        return *entry;
    }
};

} // namespace

// ==============================================================================
// Cross Section & Mean Free Path
// ==============================================================================
double RayleighScattering::computeCrossSectionPerAtom(double photonEnergy, const Material& material) const {
    if (photonEnergy <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    const int z = static_cast<int>(std::lround(material.getZ()));
    const double barn = 1.0e-24 * Units::cm * Units::cm;
    return PhotoElectricShellData::get(z).getCoherentCrossSection(photonEnergy) * barn;
}

double RayleighScattering::computeMassAttenuation(double photonEnergy, const Material& material) const {
    const double density = material.getDensity();
    if (density <= 0.0 || material.getZ() < 1.0) {
        return 0.0;
    }
    return (material.getAtomDensity() / density) * computeCrossSectionPerAtom(photonEnergy, material);
}

double RayleighScattering::getMeanFreePath(const Track& track, const Material& material) const {
    const double mu = material.getAtomDensity() * computeCrossSectionPerAtom(track.getKineticEnergy(), material);
    return mu > 0.0 ? 1.0 / mu : GeoDefs::kInfinity;
}

// ==============================================================================
// Angular Sampling (x^2 CDF + Thomson Rejection)
// ==============================================================================
double RayleighScattering::sampleCosTheta(double photonEnergy, int z) const {
    const double hc = kHcEvAngstrom * Units::eV;
    const double xMax = photonEnergy / hc;
    const double uMax = xMax * xMax;
    if (uMax <= 0.0) {
        return 1.0;
    }

    const ElementRayleighCdf& data = ElementRayleighCdf::get(z);
    const double aMax = data.evaluateCdf(uMax);
    if (aMax <= 0.0) {
        return 1.0;
    }

    double cosTheta = 1.0;
    for (int iter = 0; iter < 1000; ++iter) {
        const double u = data.sampleX2(RandomEngine::uniNumGen() * aMax, uMax);
        cosTheta = std::clamp(1.0 - 2.0 * u / uMax, -1.0, 1.0);
        if (RandomEngine::uniNumGen() <= 0.5 * (1.0 + cosTheta * cosTheta)) {
            break;
        }
    }
    return cosTheta;
}

// ==============================================================================
// PostStepDoIt: Elastic Deflection (Unchanged Energy, No Secondaries)
// ==============================================================================
void RayleighScattering::postStepDoIt(Track& track, const Material& material, Step&, std::vector<Track>&) {
    clearNumberOfInteractionLengthLeft();

    const double photonEnergy = track.getKineticEnergy();
    const int z = static_cast<int>(std::lround(material.getZ()));
    if (photonEnergy <= 0.0 || z < 1) {
        return;
    }

    const double cosTheta = sampleCosTheta(photonEnergy, z);
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();

    Vector3D newDirection(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    newDirection.rotateUz(track.getDirection());
    track.setDirection(newDirection);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] coherent scatter, cos(theta) = " << cosTheta
                  << " (E_gamma = " << photonEnergy / Units::keV << " keV)\n";
    }
}