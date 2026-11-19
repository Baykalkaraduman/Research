#pragma once

#include <string>
#include "Vector3D.h"

// ==============================================================================
// Struct: Step
// What happened during ONE step of a track. Filled by the stepping loop and read by the analysis code.
// ==============================================================================
struct Step {
    Vector3D preStepPosition;
    Vector3D postStepPosition;
    double stepLength = 0.0;
    double totalEnergyDeposit = 0.0;
    bool leftWorld = false;

    // Name of the process that limited this step ("Decay", "Ionization",
    // "Transportation", "UserLimit").
    std::string processDefinedStep;

    void addTotalEnergyDeposit(double energy) {totalEnergyDeposit += energy;}

};
