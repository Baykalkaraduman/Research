#pragma once

#include "Particle.h"
#include "Vector3D.h"

// ==============================================================================
// Class: Track
// Represents a particle propagating through the simulation geometry.
// Holds dynamic state (position, direction, time) and a pointer to static properties.
// Equivalent to G4Track in Geant4.
// ==============================================================================
enum class TrackStatus {
    trackAlive,         // Track is alive and the observed
    trackKill,          // Track is killed, stop observing
    trackStopButAlive   // Particle is at rest, waits for the next process (e.g. Annihilation)
};

class Track {
private:
    // Pointer to the immutable physical properties of the particle type.
    // 'const' ensures we cannot accidentally change the mass of an electron mid-flight!
    const Particle* particleDefinition; 

    // Dynamic Kinematic State
    Vector3D position;          // Position in 3D space (e.g., cm or mm)
    Vector3D direction;         // Direction of the particle
    double kineticEnergy;       // Kinetic energy (GeV)
    double globalTime;          // Time since the start of the event (ns)
    double trackLength = 0.0;   // Distance travelled so far by the particle (cm)
    
    // Status flag
    TrackStatus trackStatus;          // Shows whether particle is alive or not

public:
    // Constructor requires the particle definition and initial kinematics
    Track(const Particle* pDef, const Vector3D& initPos, const Vector3D& initDirection, double initKinEnergy);

    // ==============================================================================
    // Getters
    // ==============================================================================
    const Particle* getParticleDefinition() const {return particleDefinition;}
    Vector3D getPosition() const {return position;}
    Vector3D getDirection() const {return direction;}
    double getKineticEnergy() const {return kineticEnergy;}
    TrackStatus getTrackStatus() const {return trackStatus;}
    double getGlobalTime() const {return globalTime;}
    double getTrackLength() const {return trackLength;}
    
    // Tells whether the particle is alive or not
    bool isAlive() const {return trackStatus == TrackStatus::trackAlive;}

    // --------------------------------------------------------------------
    // Getters for Derived Parameters
    // --------------------------------------------------------------------
    // Calling m from particleDefinition to make writing convinient
    double getMass() const {return particleDefinition->getMass();}

    // Magnitude of Momentum, |p| = [T(T + 2m)]^(1/2)
    double getMomentumMag() const {return std::sqrt(kineticEnergy*(kineticEnergy + 2*getMass()));}

    // Momentum Vector, p^(->) = direction * |p|
    Vector3D getMomentum() const {return direction*getMomentumMag();}

    // Total Energy E = T + m
    double getTotalEnergy() const {return kineticEnergy + getMass();}

    // Beta β = |p| / E
    double getBeta() const {return getMomentumMag() / getTotalEnergy();}

    // Gamma γ = E / m
    double getGamma() const {return getTotalEnergy() / getMass();} // Don't need for massles particles and don't use it because of zero division

    // Beta * Gamma = β * γ = |p| / m
    double getBetaGamma() const {return getMomentumMag() / getMass();}

    // ==========================================================================
    // Setters / Updaters (Called by physics processes and stepping manager)
    // ==========================================================================
    void setDirection(const Vector3D& newDirection) {direction = newDirection.unit();}
    void setTrackStatus(TrackStatus status) {trackStatus = status;}
    void addTime(double deltaTime) {globalTime += deltaTime;}
    void addTrackLength(double length) {trackLength += length;}
    void deductEnergy(double energyLoss);
    void updatePosition(const Vector3D& newPos) {position = newPos;}
    void setKineticEnergy(double energy) {kineticEnergy = energy;}

    // ==========================================================================
    // Utility
    // ==========================================================================

    // Kill the track (e.g., when it hits a wall, loses all energy, or decays)
    void killTrack() {trackStatus = TrackStatus::trackKill;}

    // Display track state
    void displayState() const;
};