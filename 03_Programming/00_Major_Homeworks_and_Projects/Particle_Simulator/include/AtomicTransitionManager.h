#pragma once

#include <string>
#include <vector>

// ==============================================================================
// Class: AtomicTransitionManager
// Manages atomic subshell binding energies and radiative / Auger relaxation
// probabilities loaded from ENDF/B-VIII.0 (MF=28, MT=533) per element.
// ==============================================================================
class AtomicTransitionManager {
public:
    // --------------------------------------------------------------------------
    // Transition & Subshell Data Types
    // --------------------------------------------------------------------------
    struct RadiativeTransition {
        int originShell;
        double probability;
        double energy;     // B_s - B_o [GeV]
        double fileEnergy; // ETR from ENDF file [GeV]
    };

    struct AugerTransition {
        int originShell;
        int augerShell;
        double probability;
        double energy;     // B_s - B_o - B_a [GeV]
        double fileEnergy; // ETR from ENDF file [GeV]
    };

    struct Shell {
        int designator;
        double bindingEnergy; // [GeV]
        double electrons;
        std::vector<RadiativeTransition> radiative;
        std::vector<AugerTransition> auger;
        double totalRadiative = 0.0;
        double totalAuger = 0.0;
    };

private:
    // --------------------------------------------------------------------------
    // Member Data & Internal Construction
    // --------------------------------------------------------------------------
    int z;
    std::vector<Shell> shells;

    explicit AtomicTransitionManager(int atomicNumber);

public:
    // --------------------------------------------------------------------------
    // Element Cache & ENDF Naming Utilities
    // --------------------------------------------------------------------------
    [[nodiscard]] static const AtomicTransitionManager& get(int atomicNumber);
    [[nodiscard]] static std::string getEndfFileName(const std::string& sublibrary, int atomicNumber);
    [[nodiscard]] static std::string getShellName(int designator);

    // --------------------------------------------------------------------------
    // Subshell Accessors
    // --------------------------------------------------------------------------
    [[nodiscard]] int getZ() const {return z;}
    [[nodiscard]] int getNumberOfShells() const {return static_cast<int>(shells.size());}
    [[nodiscard]] const Shell& getShell(int index) const {return shells[index];}
    [[nodiscard]] double getBindingEnergy(int index) const {return shells[index].bindingEnergy;}
    [[nodiscard]] bool hasTransitions(int index) const {return !shells[index].radiative.empty() || !shells[index].auger.empty();}
    [[nodiscard]] int findShellIndex(int designator) const;
};