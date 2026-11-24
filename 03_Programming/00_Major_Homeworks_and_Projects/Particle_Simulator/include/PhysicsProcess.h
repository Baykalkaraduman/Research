#pragma once

#include "GeoDefs.h"
#include "Material.h"
#include "Step.h"
#include "Track.h"

#include <string>
#include <vector>

// ==============================================================================
// Abstract Base Class: PhysicsProcess
// Defines the interface for all physics interactions in the simulation.
// ==============================================================================
class PhysicsProcess {
protected:
    std::string processName;

    // What should be printed and what should not (1 for print, 0 for silent)
    int verboseLevel = 1;

public:
    // ==========================================================================
    // Constructors and Destructors
    // ==========================================================================
    explicit PhysicsProcess(std::string name) : processName(std::move(name)) {}
    virtual ~PhysicsProcess() = default;

    // ==========================================================================
    // Setup
    // ==========================================================================
    // Checks whether the process is applicable for the particle
    [[nodiscard]] virtual bool isApplicable(const Particle&) const {return true;}

    // Starts a new track when there is a new particle
    virtual void startTracking(Track& /*track*/) {}

    // Setter for verbose level
    void setVerboseLevel(int v_level) {verboseLevel = v_level;}

    // ==========================================================================
    // Getters
    // ==========================================================================
    [[nodiscard]] std::string getProcessName() const {return processName;}
    [[nodiscard]] int getVerboseLevel() const {return verboseLevel;}

    // ==========================================================================
    // ALONG-STEP (continuous) half
    // ==========================================================================
    // Max step allowed so that the continuous approximation stays valid
    [[nodiscard]] virtual double alongStepGetPhysicalInteractionLength(const Track&, const Material&) {
        return GeoDefs::kInfinity;
    }

    // Apply the continuous effect over the step that was actually taken
    virtual void alongStepDoIt(Track&, const Material&, double /*stepLength*/) {}

    // ==========================================================================
    // POST-STEP (discrete) half
    // ==========================================================================
    // Distance to the next discrete interaction of this process
    [[nodiscard]] virtual double postStepGetPhysicalInteractionLength(const Track&, const Material&) {
        return GeoDefs::kInfinity;
    }

    // Called after EVERY step, for every process:
    // the particle moves a distance stepLength, so it is that much closer to its next interaction
    virtual void subtractNumberOfInteractionLengthLeft(double /*stepLength*/) {}

    // Perform the interaction. Only called for the process that limited the step.
    virtual void postStepDoIt(Track&, const Material&, Step& /*step*/, std::vector<Track>& /*secondaries*/) {}

    // ==========================================================================
    // AT-REST half
    // ==========================================================================
    // Pure query (must NOT touch RandomEngine): can this process act at rest here?
    [[nodiscard]] virtual bool isAtRestApplicable(const Particle&, const Material&) const {return false;}

    // If true, this process owns the entire stopped-particle lifecycle in this material
    // (e.g. MuonCapture handles atomic cascade + capture/DIO competition for mu- at Z >= 1).
    [[nodiscard]] virtual bool isAtRestExclusive(const Particle&, const Material&) const {return false;}

    // Returns the sampled TIME (in ns) until the process happens, or kInfinity if inapplicable.
    [[nodiscard]] virtual double atRestGetPhysicalInteractionLength(const Track&, const Material&) {
        return GeoDefs::kInfinity;
    }

    // What happens at rest. Must kill the track.
    virtual void atRestDoIt(Track&, const Material&, Step& /*step*/, std::vector<Track>& /*secondaries*/) {}
};