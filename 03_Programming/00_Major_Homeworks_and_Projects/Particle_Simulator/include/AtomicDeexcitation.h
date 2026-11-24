#pragma once

#include "Track.h"
#include "Vector3D.h"

#include <vector>

class Particle;

// ==============================================================================
// Class: AtomicDeexcitation
// Simulates atomic relaxation cascades (radiative X-ray fluorescence and
// non-radiative Auger emission) following an inner-shell vacancy, returning
// unemitted binding energy as local energy deposit.
// ==============================================================================
class AtomicDeexcitation {
private:
    // --------------------------------------------------------------------------
    // Member Data & Internal Helpers
    // --------------------------------------------------------------------------
    const Particle* gammaDef;
    const Particle* electronDef;
    bool fluorescence = true;
    bool auger = false;

    [[nodiscard]] static Vector3D sampleIsotropicDirection();

public:
    // --------------------------------------------------------------------------
    // Construction & Configuration
    // --------------------------------------------------------------------------
    AtomicDeexcitation(const Particle* gamma, const Particle* electron) : gammaDef(gamma), electronDef(electron) {}

    void setFluorescence(bool on) {fluorescence = on;}
    void setAuger(bool on) {auger = on;}

    [[nodiscard]] bool isFluorescenceActive() const {return fluorescence;}
    [[nodiscard]] bool isAugerActive() const {return auger;}

    // --------------------------------------------------------------------------
    // Cascade Generation
    // --------------------------------------------------------------------------
    // Relaxes a vacancy in 'shellIndex' of element Z at 'position', appends
    // emitted secondaries, and returns the local energy deposit.
    [[nodiscard]] double generateParticles(int z, int shellIndex, const Vector3D& position, std::vector<Track>& secondaries) const;
};