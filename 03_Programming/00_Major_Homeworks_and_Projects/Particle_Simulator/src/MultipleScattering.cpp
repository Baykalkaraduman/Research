#include "MultipleScattering.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"

#include <algorithm> // std::max
#include <cmath>
#include <iostream>

// ==============================================================================
// Highland formula, squared
// ==============================================================================
double MultipleScattering::computeTheta0Squared(const Track& track, const Material& material, double length) const {
    const double radLength = material.getRadLength();
    if (length <= 0.0 || radLength >= GeoDefs::kInfinity) {
        return 0.0;
    }
    const double charge = track.getParticleDefinition()->getCharge();
    const double beta = track.getBeta();
    const double momentum = track.getMomentumMag();

    const double xOverX0 = length / radLength;
    const double correction = 1.0 + 0.038 * std::log(xOverX0 * charge * charge / (beta * beta));
    const double theta0 = PhysicalConstants::highlandConstant / (beta * momentum) * std::abs(charge)
                        * std::sqrt(xOverX0) * correction;
    return theta0 * theta0;
}

double MultipleScattering::computeTheta0(const Track& track, const Material& material, double thickness) const {
    return std::sqrt(computeTheta0Squared(track, material, thickness));
}

// ==============================================================================
// Along-step DoIt: deflect the direction
// ==============================================================================
void MultipleScattering::alongStepDoIt(Track& track, const Material& material, double stepLength) {
    if (material.getRadLength() >= GeoDefs::kInfinity) {
        return; // vacuum: nothing to scatter on
    }

    // New material (or first step of the track): a new layer starts
    if (&material != currentMaterial) {
        currentMaterial = &material;
        pathLengthInMaterial = 0.0;
    }

    // Variance added by THIS step (see the header for why it is a difference).
    // max(0, ...): for absurdly thin layers (x/X0 < 1e-11) the log term passes -1
    // and theta0^2(L) is no longer increasing; never happens, but must not give NaN.
    const double before = computeTheta0Squared(track, material, pathLengthInMaterial);
    pathLengthInMaterial += stepLength;
    const double after = computeTheta0Squared(track, material, pathLengthInMaterial);
    const double theta0 = std::sqrt(std::max(after - before, 0.0));

    // Two projected angles, independent Gaussians of width theta0 (PDG eq. 34.16/34.17).
    // The space angle is theta = sqrt(thetaX^2 + thetaY^2), the azimuth atan2.
    const double thetaX = RandomEngine::randGauss(0.0, theta0);
    const double thetaY = RandomEngine::randGauss(0.0, theta0);
    const double theta = std::sqrt(thetaX * thetaX + thetaY * thetaY);
    const double phi = std::atan2(thetaY, thetaX);

    // New direction in the frame where the old direction is the z axis, then
    // rotated to the global frame.
    Vector3D newDirection(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
    newDirection.rotateUz(track.getDirection());
    track.setDirection(newDirection);

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] theta0 = " << theta0 * 1e3 << " mrad, deflected by "
                  << theta * 1e3 << " mrad\n";
    }
}
