#pragma once

#include "DecayChannel.h"

#include <memory>
#include <vector>

// ==============================================================================
// Class: DecayTable
// All decay channels of one particle type, with their branching ratios.
// A Particle points to its table, the decay process asks it for a channel.
// ==============================================================================
class DecayTable {
private:
    std::vector<std::unique_ptr<DecayChannel>> channels;

public:
    void insert(std::unique_ptr<DecayChannel> channel) {channels.push_back(std::move(channel));}
    int entries() const {return static_cast<int>(channels.size());}
    const DecayChannel* getDecayChannel(int index) const {return channels[index].get();}

    // Sum of the branching ratios (should be 1; the selection divides by it anyway)
    double getSumOfBR() const;

    // Pick a channel at random with probability BR_i / sum(BR). -> DecayTable.cpp
    // (G4DecayTable::SelectADecayChannel)
    const DecayChannel* selectADecayChannel() const;
};
