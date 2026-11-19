#pragma once

#include "Units.h"

namespace GeoDefs {

    // Predefined large number for physical processes with unlimited step length
    // e.g. stable particles that missed the target etc.
    constexpr double kInfinity = 9.9e+99;

    // Moves particle inside/outside the material volume so it cannot stuck on the surface.
    constexpr double kCarTolerance = 1.0e-10 * Units::cm;
}