#pragma once

#include "Volume.h"
#include "Track.h"
#include "Step.h"
#include "AnalysisManager.h"
#include "PhysicsProcess.h"

#include <vector>
#include <memory>
#include <functional>

// ==============================================================================
// Class: RunManager
// Manages geometry, physics processes, and the main event stepping loop.
// ==============================================================================
class RunManager {
private:
    std::unique_ptr<Volume> worldVolume;
    std::vector<std::unique_ptr<Volume>> geometryDB; // daughters: inside the world, must not overlap
    std::vector<std::unique_ptr<PhysicsProcess>> physicsList;
    AnalysisManager* analysisManager = nullptr;
    int verboseLevel = 1;
    std::function<void(const Track&, const Step&)> userSteppingAction;
    double energyCut = 0.0;
    int numberOfTracks = 0;
    std::size_t maxStackSize = 0;

    // ==========================================================================
    // Internal stepping machinery
    // ==========================================================================
    void trackParticle(Track& track, std::vector<Track>& trackStack, double maxStepLength, int maxSteps);
    Step processOneStep(Track& track, std::vector<Track>& secondaries, double maxStepLength);
    Step processAtRest(Track& track, std::vector<Track>& secondaries);
    Step applyEnergyCut(Track& track);

    bool isBelowEnergyCut(const Track& track) const;
    bool hasAtRestProcess(const Particle& particle);

    // Distance to the next geometry boundary along the track direction
    double computeGeometryStep(const Vector3D& point, const Vector3D& direction) const;

public:
    RunManager() = default;
    ~RunManager() = default;

    // ==========================================================================
    // Setting and Updating Methods
    // ==========================================================================
    void setWorldVolume(std::unique_ptr<Volume> world) {worldVolume = std::move(world);}
    void addVolume(std::unique_ptr<Volume> vol);
    void addPhysicsProcess(std::unique_ptr<PhysicsProcess> process);
    void setAnalysisManager(AnalysisManager* am) {analysisManager = am;}
    void setVerboseLevel(int level);
    void setUserSteppingAction(std::function<void(const Track&, const Step&)> action) {userSteppingAction = std::move(action);}
    void setEnergyCut(double cut) {energyCut = cut;}
    double getEnergyCut() const {return energyCut;}
    int getNumberOfTracks() const {return numberOfTracks;}
    std::size_t getMaxStackSize() const {return maxStackSize;}

    // Daughter volume containing the point, else the world, else nullptr (outside the world)
    Volume* getVolumeAt(const Vector3D& point) const;

    // ==========================================================================
    // Simulation Execution
    // ==========================================================================
    void simulateEvent(const Track& primaryTrack, double maxStepLength, int maxSteps);
};
