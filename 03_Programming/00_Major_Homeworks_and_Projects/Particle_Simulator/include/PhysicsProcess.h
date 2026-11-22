#pragma once

#include "Track.h"
#include "Material.h"
#include "GeoDefs.h"
#include "Step.h"

#include <string>
#include <vector>


// ==============================================================================
// Abstract Base Class: PhysicsProcess
// Defines the interface for all physics interactions in the simulation.
// ==============================================================================
class PhysicsProcess {
protected:
    std::string processName;

    // What should be printed and what should not
    // 1 for print, 0 for silent
    int verboseLevel = 1;

public:
    // ==========================================================================
    // Constructors and Destructors
    // ==========================================================================

    // Constructor
    explicit PhysicsProcess(std::string name) : processName(std::move(name)) {}

    // Virtual Destructor
    virtual ~PhysicsProcess() = default;

    // ==========================================================================
    // Setup
    // ==========================================================================
    // Checks whether the process is applicable for the particle
    virtual bool isApplicable(const Particle&) const {return true;}

    // Starts a new track when there is a new particle
    virtual void startTracking(Track& /*track*/) {}

    // Setter for verbose level
    void setVerboseLevel(int v_level) {verboseLevel = v_level;}
    
    // ==========================================================================
    // Getters
    // ==========================================================================

    std::string getProcessName() const {return processName;}
    int getVerboseLevel() const {return verboseLevel;}

    // ==========================================================================
    // ALONG-STEP (continuous) half
    // ==========================================================================
    // Max step allowed so that the continuous approximation stays valid
    virtual double alongStepGetPhysicalInteractionLength(const Track&, const Material&) {return GeoDefs::kInfinity;}

    // Apply the continuous effect over the step that was actually taken
    virtual void alongStepDoIt(Track&, const Material&, double /*stepLength*/) {}

    // ==========================================================================
    // POST-STEP (discrete) half
    // ==========================================================================
    // Distance to the next discrete interaction of this process
    virtual double postStepGetPhysicalInteractionLength(const Track&, const Material&) {return GeoDefs::kInfinity;}

    // Called after EVERY step, for every process
    // The particle moves a distance`stepLength`, so it is that much closer to its next interaction
    virtual void subtractNumberOfInteractionLengthLeft(double /*stepLength*/) {}

    // Perform the interaction. Only called for the process that limited the step.
    virtual void postStepDoIt(Track&, const Material&, Step& /*step*/, std::vector<Track>& /*secondaries*/) {}

    // ==========================================================================
    // AT-REST half
    // ==========================================================================
    // Only for a track with status "trackStopButAlive". Returns the mean TIME until the process happens
    virtual double atRestGetPhysicalInteractionLength(const Track&) {return GeoDefs::kInfinity;}

    // What happens at rest. Must kill the track.
    virtual void atRestDoIt(Track&, const Material&, Step& /*step*/, std::vector<Track>& /*secondaries*/) {}

};
