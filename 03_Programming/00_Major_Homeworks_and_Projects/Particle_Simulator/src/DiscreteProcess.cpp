#include "DiscreteProcess.h"
#include "RandomEngine.h"

#include <cmath>

// ==============================================================================
// GPIL: distance to the next interaction of this process
// ==============================================================================
double DiscreteProcess::postStepGetPhysicalInteractionLength(const Track& track, const Material& material) {
    // 1. New track, or the process just interacted: sample how many mean free
    //    paths the particle travels before the next interaction.
    //    flat() is in [0, 1), so 1 - flat() is in (0, 1] and log() never sees 0.
    if (numberOfInteractionLengthLeft < 0.0) {
        numberOfInteractionLengthLeft = -std::log(1.0 - RandomEngine::uniNumGen());
    }

    // 2. Mean free path HERE (depends on energy and material at the pre-step point)
    currentInteractionLength = getMeanFreePath(track, material);
    if (currentInteractionLength >= GeoDefs::kInfinity) {
        return GeoDefs::kInfinity;
    }

    // 3. Convert "mean free paths left" into a distance
    return numberOfInteractionLengthLeft * currentInteractionLength;
}

// ==============================================================================
// Bookkeeping after every step (also steps this process did not limit)
// ==============================================================================
void DiscreteProcess::subtractNumberOfInteractionLengthLeft(double stepLength) {
    if (currentInteractionLength >= GeoDefs::kInfinity) {
        return; // the process cannot happen here, so nothing is "used up"
    }
    numberOfInteractionLengthLeft -= stepLength / currentInteractionLength;

    // Rounding can push it slightly below 0. 0 means "interact at the next step".
    // (Must not become negative: negative means "sample again".)
    if (numberOfInteractionLengthLeft < 0.0) {
        numberOfInteractionLengthLeft = 0.0;
    }
}
