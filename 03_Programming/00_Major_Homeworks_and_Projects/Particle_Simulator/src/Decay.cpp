#include "Decay.h"
#include "DecayTable.h"
#include "RandomEngine.h"

#include <algorithm>
#include <iostream>
#include <cmath>

// ==============================================================================
// The decay time, drawn once per track
// ==============================================================================
// P(tau_decay > t) = exp(-t / tau)  ->  tau_decay = -tau ln(u),  u uniform in (0, 1].
// uniNumGen() is in [0, 1), so 1 - uniNumGen() is in (0, 1] and log() never sees 0.
// A value given from outside (an event generator, a test) is kept, like Geant4
// keeps a pre-assigned decay time.
void Decay::startTracking(Track& track) {
    const Particle& particle = *track.getParticleDefinition();
    if (!isApplicable(particle) || track.getPreAsgnDecayPrpT() >= 0.0) {
        return;
    }
    track.setPreAsgnDecayPrpT(-particle.getLifetime() * std::log(1.0 - RandomEngine::uniNumGen()));
}

double Decay::getRemainingProperTime(const Track& track) {
    const double preAssigned = track.getPreAsgnDecayPrpT();
    if (preAssigned < 0.0) {
        // Not drawn (RunManager::hasAtRestProcess asks with a fresh track): the
        // mean lifetime, only to say "yes, this process acts at rest"
        return track.getParticleDefinition()->getLifetime();
    }
    // Rounding in the proper-time sum can overshoot by ~1e-16 relative: clamp at 0
    return std::max(0.0, preAssigned - track.getProperTime());
}

// ==============================================================================
// The decay itself (in flight or at rest)
// ==============================================================================
void Decay::decayIt(Track& track, std::vector<Track>& secondaries) {
    const Particle& particle = *track.getParticleDefinition();
    const DecayChannel* channel = particle.getDecayTable()->selectADecayChannel();

    if (verboseLevel > 0) {
        std::cout << "    [" << processName << "] *** DECAY TRIGGERED! *** " << channel->getDescription()
                  << "  (proper time " << track.getProperTime() << " ns, "
                  << (track.getKineticEnergy() > 0.0 ? "in flight" : "at rest") << ")\n";
    }

    std::vector<LorentzVector> products = channel->decayIt(particle.getMass());

    // Boost from the parent's rest frame to the lab. At rest the boost is the
    // identity (beta = 0), LorentzVector::boost returns at once.
    const double parentEnergy = track.getTotalEnergy();
    const Vector3D beta = track.getMomentum() / parentEnergy;
    const double gamma = parentEnergy / particle.getMass();

    for (int i = 0; i < channel->getNumberOfDaughters(); ++i) {
        LorentzVector& p4 = products[i];
        p4.boost(beta, gamma);
        const Particle* daughter = channel->getDaughter(i);
        // E - m can come out as -1e-17 for a daughter almost at rest: clamp
        const double kineticEnergy = std::max(0.0, p4.E() - daughter->getMass());
        secondaries.emplace_back(daughter, track.getPosition(), p4.vect(), kineticEnergy);
        // [STAGE 3] The parent's time is given to EVERY secondary by RunManager.

        if (verboseLevel > 0) {
            std::cout << "    [" << processName << "]   -> " << daughter->getName()
                      << "  T = " << kineticEnergy / Units::MeV << " MeV\n";
        }
    }

    // ==========================================================================
    // TERMINATE PARENT TRACK
    // ==========================================================================
    track.killTrack();
}