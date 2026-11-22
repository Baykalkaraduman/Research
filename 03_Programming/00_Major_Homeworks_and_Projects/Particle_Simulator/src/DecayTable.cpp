#include "DecayTable.h"
#include "RandomEngine.h"

// ==============================================================================
// Decay Table
// ==============================================================================
double DecayTable::getSumOfBR() const {
    double sum = 0.0;
    for (const auto& channel : channels) {
        sum += channel->getBR();
    }
    return sum;
}

// One uniform number u in [0, sum(BR)); walk through the channels adding up the BRs until the running sum passes u (inverse of the cumulative distribution).
const DecayChannel* DecayTable::selectADecayChannel() const {
    if (channels.empty()) {
        return nullptr;
    }
    const double u = RandomEngine::uniNumGen() * getSumOfBR();
    double running = 0.0;
    for (const auto& channel : channels) {
        running += channel->getBR();
        if (u < running) {
            return channel.get();
        }
    }
    return channels.back().get(); // u == sum by rounding
}
