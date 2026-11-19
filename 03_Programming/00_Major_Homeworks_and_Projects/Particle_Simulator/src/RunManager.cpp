#include "RunManager.h"
#include "PhysicalConstants.h" // Required for speed of light
#include "GeoDefs.h"

#include <iostream>
#include <algorithm>  // std::min
#include <stdexcept>  // std::runtime_error

// ==========================================================================
// Geometry & Physics Registration
// ==========================================================================

void RunManager::addVolume(std::unique_ptr<Volume> vol) {
    geometryDB.push_back(std::move(vol));
}

void RunManager::addPhysicsProcess(std::unique_ptr<PhysicsProcess> process) {
    physicsList.push_back(std::move(process));
}

void RunManager::setVerboseLevel(int level) {
    verboseLevel = level;
    for (const auto& process : physicsList) {
        process->setVerboseLevel(level);
    }
}

Volume* RunManager::getVolumeAt(const Vector3D& point) const {
    for (const auto& vol : geometryDB) {
        if (vol->isInside(point)) {
            return vol.get();
        }
    }
    // Not in any daughter -> in the world (if still inside it)
    if (worldVolume && worldVolume->isInside(point)) {
        return worldVolume.get();
    }
    return nullptr; // outside the world
}

// ==========================================================================
// Transportation: distance to the next boundary
// ==========================================================================
double RunManager::computeGeometryStep(const Vector3D& point, const Vector3D& direction) const {
    // Inside a daughter: the only boundary ahead is its own surface
    for (const auto& vol : geometryDB) {
        if (vol->isInside(point)) {
            return vol->distanceToOut(point, direction);
        }
    }

    // In the world, outside all daughters: either we enter the nearest daughter
    // on our way, or we reach the world's own wall
    double distance = worldVolume->distanceToOut(point, direction);
    for (const auto& vol : geometryDB) {
        distance = std::min(distance, vol->distanceToIn(point, direction));
    }
    return distance;
}

// ==========================================================================
// Main Event Simulation Loop (Track Stack & Secondaries)
// ==========================================================================

void RunManager::simulateEvent(const Track& primaryTrack, double maxStepLength, int maxSteps) {
    if (!worldVolume) {
        throw std::runtime_error("RunManager: no world volume set (call setWorldVolume first)");
    }

    if (verboseLevel > 0) std::cout << "\n=== RUN MANAGER: EVENT START ===\n";

    // The Track Stack (LIFO buffer)
    std::vector<Track> trackStack;
    trackStack.push_back(primaryTrack);

    int trackCounter = 1;
    numberOfTracks = 0;
    maxStackSize = 1;

    // Main Event Loop: Runs until no particles are left in the stack
    while (!trackStack.empty()) {
        // Pop the last added track from the stack
        Track currentTrack = trackStack.back();
        trackStack.pop_back();

        if (verboseLevel > 0) {
            std::cout << "\n==================================================\n";
            std::cout << ">>> TRACKING PARTICLE " << trackCounter
                      << " (" << currentTrack.getParticleDefinition()->getName() << ") >>>\n";
            std::cout << "==================================================\n";
        }
        trackCounter++;
        numberOfTracks++;

        trackParticle(currentTrack, trackStack, maxStepLength, maxSteps);
        maxStackSize = std::max(maxStackSize, trackStack.size());
    }
    if (verboseLevel > 0) std::cout << "\n=== RUN MANAGER: EVENT COMPLETED (No tracks left) ===\n";
}

// ==========================================================================
// One track, step by step
// ==========================================================================
void RunManager::trackParticle(Track& track, std::vector<Track>& trackStack, double maxStepLength, int maxSteps) {
    if (getVolumeAt(track.getPosition()) == nullptr) {
        if (verboseLevel > 0) std::cout << "    [Warning] Track starts outside the world. Killed.\n";
        track.killTrack();
        return;
    }

    // Every discrete process forgets what it sampled for the previous track
    for (const auto& process : physicsList) {
        process->startTracking();
    }

    int stepCounter = 0;
    while (track.getTrackStatus() != TrackStatus::trackKill && stepCounter < maxSteps) {
        if (verboseLevel > 0) {
            std::cout << "\n--- Step " << stepCounter << " ---\n";
            track.displayState();

            Volume* currentVolume = getVolumeAt(track.getPosition());
            std::cout << "    [Geometry] Particle is inside: " << currentVolume->getName()
                      << " (Material: " << currentVolume->getMaterial()->getName() << ")\n";
        }
        // at rest -> at-rest process;  below the energy cut -> cut;  else -> normal step
        std::vector<Track> stepSecondaries;
        Step step;
        if (track.getTrackStatus() == TrackStatus::trackStopButAlive) {
            step = processAtRest(track, stepSecondaries);
        } else if (isBelowEnergyCut(track)) {
            step = applyEnergyCut(track);
        } else {
            step = processOneStep(track, stepSecondaries, maxStepLength);
        }

        // Stopped, but nothing acts at rest on this particle (e-, mu-, ...):
        // it simply stays where it is.
        if (track.getTrackStatus() == TrackStatus::trackStopButAlive && !hasAtRestProcess(*track.getParticleDefinition())) {
            track.killTrack();
        }

        if (verboseLevel > 0) {
            std::cout << "    [Step] Length: " << step.stepLength << " cm | Defined by: "
                      << step.processDefinedStep << "\n";
        }

        // User hook
        if (userSteppingAction) {
            userSteppingAction(track, step);
        }

        // =========================================================
        // ANALYSIS
        // =========================================================
        if (analysisManager != nullptr) {
            if (step.totalEnergyDeposit > 0.0) {
                analysisManager->fillEnergyDeposition(step.totalEnergyDeposit);
            }
            if (!track.isAlive() && step.processDefinedStep == "Decay") {
                analysisManager->fillDecayPosition(track.getPosition().z);
                analysisManager->fillTimeOfFlight(track.getGlobalTime());
            }
        }

        // Push secondaries to the track stack
        for (auto& sec : stepSecondaries) {
            sec.addTime(track.getGlobalTime());
            trackStack.push_back(sec);
        }
        stepCounter++;
    }

    if (verboseLevel > 0) {
        std::cout << "\n  [Final State]\n";
        track.displayState();
    }
    if (track.isAlive()) {
        // Printed even when silent: it means the result is incomplete
        std::cout << "    [Warning] maxSteps reached, track abandoned while still alive.\n";
    }
}

// ==========================================================================
// One step
// ==========================================================================
//   1. GPIL:   every process proposes a step length, geometry too -> take the shortest
//   2. Move:   position and time
//   3. AlongStepDoIt of ALL processes (continuous effects over the step)
//   4. Bookkeeping: every discrete process "uses up" stepLength / lambda
//   5. PostStepDoIt of the ONE discrete process that won (if it won)
//   6. Left the world? -> kill
Step RunManager::processOneStep(Track& track, std::vector<Track>& secondaries, double maxStepLength) {
    Step step;
    step.preStepPosition = track.getPosition();

    // Material at the PRE-step point is used for the whole step. Because a step
    // never crosses a boundary (step 1), the whole step really is in this material.
    const Material& material = *getVolumeAt(track.getPosition())->getMaterial();
    const Particle& particle = *track.getParticleDefinition();

    // =========================================================
    // 1. GPIL: WHO limits the step?
    // =========================================================
    double stepLength = maxStepLength;        // user limit
    step.processDefinedStep = "UserLimit";
    PhysicsProcess* postStepProcess = nullptr; // discrete process that will act at the end

    // Discrete processes: distance to their next interaction
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        double proposed = process->postStepGetPhysicalInteractionLength(track, material);
        if (proposed < stepLength) {
            stepLength = proposed;
            postStepProcess = process.get();
            step.processDefinedStep = process->getProcessName();
        }
    }

    // Continuous processes: max step for their approximation to stay valid
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        double proposed = process->alongStepGetPhysicalInteractionLength(track, material);
        if (proposed < stepLength) {
            stepLength = proposed;
            postStepProcess = nullptr; // a continuous limit: no discrete interaction this step
            step.processDefinedStep = process->getProcessName();
        }
    }

    // Geometry: never cross a boundary inside a step
    double geometryStep = computeGeometryStep(track.getPosition(), track.getDirection());
    if (geometryStep <= stepLength) {
        // Stop ON the boundary, then kCarTolerance across it (see GeoDefs.h)
        stepLength = geometryStep + GeoDefs::kCarTolerance;
        postStepProcess = nullptr;
        step.processDefinedStep = "Transportation";
    }

    // =========================================================
    // 2. MOVE: position and time
    // =========================================================
    // Uses the pre-step velocity. Fine while steps lose little energy, which is
    // exactly what Ionization's step limit guarantees.
    double dt = stepLength / (track.getBeta() * PhysicalConstants::c);
    track.updatePosition(track.getPosition() + track.getDirection() * stepLength);
    track.addTime(dt);
    track.addTrackLength(stepLength);

    // =========================================================
    // 3. ALONG-STEP DoIt: continuous processes (energy loss, ...)
    // =========================================================
    double energyBefore = track.getKineticEnergy();
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        process->alongStepDoIt(track, material, stepLength);
        if (!track.isAlive()) {
            // Stopped by energy loss: THIS process ended the track
            step.processDefinedStep = process->getProcessName();
            break;
        }
    }
    step.totalEnergyDeposit = energyBefore - track.getKineticEnergy();

    // =========================================================
    // 4. Bookkeeping: all discrete processes move closer to their interaction
    // =========================================================
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        process->subtractNumberOfInteractionLengthLeft(stepLength);
    }

    // =========================================================
    // 5. POST-STEP DoIt: only the discrete process that won, if any
    // =========================================================
    if (postStepProcess != nullptr && track.isAlive()) {
        postStepProcess->postStepDoIt(track, material, step, secondaries);
    }

    // =========================================================
    // 6. Left the world?
    // =========================================================
    if (track.isAlive() && getVolumeAt(track.getPosition()) == nullptr) {
        if (verboseLevel > 0) std::cout << "    [Transportation] Particle left the world. Killed.\n";
        track.killTrack();
        step.leftWorld = true;  // Kinetic Energy leaves the world as well
    }

    step.stepLength = stepLength;
    step.postStepPosition = track.getPosition();
    return step;
}
// ==========================================================================
// At rest
// ==========================================================================
// Among the processes that act at rest, the one with the shortest mean time
// wins (for us: only e+ annihilation, time 0). A zero-length step.
Step RunManager::processAtRest(Track& track, std::vector<Track>& secondaries) {
    Step step;
    step.preStepPosition = track.getPosition();
    step.postStepPosition = track.getPosition();
    step.processDefinedStep = "AtRest";

    const Material& material = *getVolumeAt(track.getPosition())->getMaterial();
    const Particle& particle = *track.getParticleDefinition();

    PhysicsProcess* atRestProcess = nullptr;
    double shortestTime = GeoDefs::kInfinity;
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        const double time = process->atRestGetPhysicalInteractionLength(track);
        if (time < shortestTime) {
            shortestTime = time;
            atRestProcess = process.get();
        }
    }

    if (atRestProcess != nullptr) {
        track.addTime(shortestTime);
        atRestProcess->atRestDoIt(track, material, step, secondaries);
        step.processDefinedStep = atRestProcess->getProcessName();
    }
    track.killTrack(); // an at-rest process always ends the track
    return step;
}

bool RunManager::hasAtRestProcess(const Particle& particle) {
    for (const auto& process : physicsList) {
        if (!process->isApplicable(particle)) continue;
        if (process->atRestGetPhysicalInteractionLength(Track(&particle, Vector3D(), Vector3D(0, 0, 1), 0.0)) < GeoDefs::kInfinity) {
            return true;
        }
    }
    return false;
}

// ==========================================================================
// Energy cut
// ==========================================================================
// Only e-, e+ and gamma: the particles of the shower, of which there are
// thousands. (Muons stop by their own energy loss, see Ionization.)
bool RunManager::isBelowEnergyCut(const Track& track) const {
    if (energyCut <= 0.0 || track.getKineticEnergy() >= energyCut) {
        return false;
    }
    const int pdg = track.getParticleDefinition()->getPdgId();
    return pdg == 11 || pdg == -11 || pdg == 22;
}

// The kinetic energy stays HERE. The particle stops: e- and gamma are then
// killed (no at-rest process), a positron still annihilates at rest, so its
// two 511 keV photons are not lost.
Step RunManager::applyEnergyCut(Track& track) {
    Step step;
    step.preStepPosition = track.getPosition();
    step.postStepPosition = track.getPosition();
    step.processDefinedStep = "EnergyCut";
    step.addTotalEnergyDeposit(track.getKineticEnergy());
    track.deductEnergy(track.getKineticEnergy()); // -> T = 0, trackKill
    return step;
}

