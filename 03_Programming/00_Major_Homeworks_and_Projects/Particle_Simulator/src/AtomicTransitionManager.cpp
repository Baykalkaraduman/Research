#include "AtomicTransitionManager.h"
#include "EndfSection.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

// ==============================================================================
// File & Subshell Naming Utilities
// ==============================================================================
std::string AtomicTransitionManager::getEndfFileName(const std::string& sublibrary, int atomicNumber) {
    static constexpr const char* symbols[] = {
        "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"};
    if (atomicNumber < 1 || atomicNumber > 100) {
        throw std::runtime_error("AtomicTransitionManager: no element Z = " + std::to_string(atomicNumber));
    }

    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "%s/%s-%03d_%s_000.endf",
                  TOYSIM_DATA_DIR, sublibrary.c_str(), atomicNumber, symbols[atomicNumber]);
    return buffer;
}

std::string AtomicTransitionManager::getShellName(int designator) {
    static constexpr struct { char letter; int count; } groups[] = {
        {'K', 1}, {'L', 3}, {'M', 5}, {'N', 7}, {'O', 9}, {'P', 11}, {'Q', 13}};
    int first = 1;
    for (const auto& g : groups) {
        if (designator < first + g.count) {
            return g.count == 1 ? std::string(1, g.letter) : g.letter + std::to_string(designator - first + 1);
        }
        first += g.count;
    }
    return "#" + std::to_string(designator);
}

int AtomicTransitionManager::findShellIndex(int designator) const {
    for (std::size_t i = 0; i < shells.size(); ++i) {
        if (shells[i].designator == designator) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// ==============================================================================
// Construction (ENDF/B-VIII.0 Atomic Relaxation, MF=28 MT=533)
// ==============================================================================
AtomicTransitionManager::AtomicTransitionManager(int atomicNumber) : z(atomicNumber) {
    EndfSection section(getEndfFileName("atom", z), 28, 533);
    const EndfSection::Cont head = section.readCont();
    const auto nShells = static_cast<std::size_t>(std::max(0L, head.n1));

    struct RawTransition {
        int subj;
        int subk;
        double fileEnergy;
        double probability;
    };

    shells.reserve(nShells);
    std::vector<std::vector<RawTransition>> raw(nShells);

    for (std::size_t s = 0; s < nShells; ++s) {
        const EndfSection::List list = section.readList();
        const auto nTransitions = static_cast<std::size_t>(std::max(0L, list.head.n2));
        if (list.values.size() < 6 * (nTransitions + 1)) {
            throw std::runtime_error("AtomicTransitionManager: truncated shell record for Z = " + std::to_string(z));
        }

        shells.push_back({static_cast<int>(std::lround(list.head.c1)),
                          list.values[0] * Units::eV,
                          list.values[1],
                          {},
                          {}});

        auto& transitions = raw[s];
        transitions.reserve(nTransitions);
        for (std::size_t t = 0; t < nTransitions; ++t) {
            const std::size_t k = 6 * (t + 1);
            transitions.push_back({static_cast<int>(std::lround(list.values[k])),
                                   static_cast<int>(std::lround(list.values[k + 1])),
                                   list.values[k + 2] * Units::eV,
                                   list.values[k + 3]});
        }
    }

    // Recompute transition energies from this file's binding energies (EPICS2017)
    // rather than using ETR (computed from older EADL binding energies) so that
    // energy is conserved exactly during relaxation cascades. Energetically
    // forbidden transitions (energy <= 0) are skipped and excluded from the totals.
    for (std::size_t s = 0; s < shells.size(); ++s) {
        const double bs = shells[s].bindingEnergy;
        std::size_t nRad = 0;
        std::size_t nAug = 0;
        for (const RawTransition& t : raw[s]) {
            if (t.subk == 0) ++nRad;
            else ++nAug;
        }
        shells[s].radiative.reserve(nRad);
        shells[s].auger.reserve(nAug);

        for (const RawTransition& t : raw[s]) {
            const int origin = findShellIndex(t.subj);
            if (origin < 0) {
                continue;
            }
            const double bo = shells[origin].bindingEnergy;
            if (t.subk == 0) {
                const double energy = bs - bo;
                if (energy <= 0.0) {
                    continue;
                }
                shells[s].radiative.push_back({origin, t.probability, energy, t.fileEnergy});
                shells[s].totalRadiative += t.probability;
            } else {
                const int augerShell = findShellIndex(t.subk);
                if (augerShell < 0) {
                    continue;
                }
                const double energy = bs - bo - shells[augerShell].bindingEnergy;
                if (energy <= 0.0) {
                    continue;
                }
                shells[s].auger.push_back({origin, augerShell, t.probability, energy, t.fileEnergy});
                shells[s].totalAuger += t.probability;
            }
        }
    }
}

// ==============================================================================
// Element Cache
// ==============================================================================
const AtomicTransitionManager& AtomicTransitionManager::get(int atomicNumber) {
    if (atomicNumber < 1 || atomicNumber > 100) {
        throw std::runtime_error("AtomicTransitionManager: no element Z = " + std::to_string(atomicNumber));
    }

    static std::array<std::unique_ptr<AtomicTransitionManager>, 101> managers;
    auto& entry = managers[static_cast<std::size_t>(atomicNumber)];
    if (!entry) {
        entry.reset(new AtomicTransitionManager(atomicNumber));
    }
    return *entry;
}