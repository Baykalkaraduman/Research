#pragma once

#include "Vector3D.h"

#include <string>
#include <vector>

class Track; // Forward declaration for the class

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

    // The Energy that is left in the interaction point
    void addTotalEnergyDeposit(double energy) {totalEnergyDeposit += energy;}

    // The step that generates secondaries
    const std::vector<Track>* secondaryInCurrentStep = nullptr;

};
