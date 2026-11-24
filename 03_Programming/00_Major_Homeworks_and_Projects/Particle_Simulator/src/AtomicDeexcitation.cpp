#include "AtomicDeexcitation.h"
#include "AtomicTransitionManager.h"
#include "RandomEngine.h"

#include <algorithm>
#include <array>
#include <cmath>

// ==============================================================================
// Internal Helpers
// ==============================================================================
Vector3D AtomicDeexcitation::sampleIsotropicDirection() {
    const double cosTheta = 2.0 * RandomEngine::uniNumGen() - 1.0;
    const double sinTheta = std::sqrt(std::max(0.0, (1.0 - cosTheta) * (1.0 + cosTheta)));
    const double phi = 2.0 * M_PI * RandomEngine::uniNumGen();
    return Vector3D(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
}

// ==============================================================================
// Vacancy Cascade
// ==============================================================================
double AtomicDeexcitation::generateParticles(int z, int shellIndex, const Vector3D& position,
                                             std::vector<Track>& secondaries) const {
    const AtomicTransitionManager& atom = AtomicTransitionManager::get(z);
    if (shellIndex < 0 || shellIndex >= atom.getNumberOfShells()) {
        return 0.0;
    }

    const double initialBinding = atom.getBindingEnergy(shellIndex);
    if (!fluorescence) {
        return initialBinding;
    }

    double emitted = 0.0;
    std::array<int, 128> vacancies{};
    int top = 0;
    vacancies[top++] = shellIndex;

    while (top > 0) {
        const int s = vacancies[--top];
        if (!atom.hasTransitions(s)) {
            continue;
        }

        const AtomicTransitionManager::Shell& shell = atom.getShell(s);
        const double totalProb = shell.totalRadiative + shell.totalAuger;
        if (totalProb <= 0.0) {
            continue;
        }

        double r = RandomEngine::uniNumGen() * totalProb;

        // 1. Radiative transitions (X-ray fluorescence)
        const AtomicTransitionManager::RadiativeTransition* radChosen = nullptr;
        for (const auto& t : shell.radiative) {
            r -= t.probability;
            if (r <= 0.0) {
                radChosen = &t;
                break;
            }
        }
        if (radChosen == nullptr && shell.auger.empty() && !shell.radiative.empty()) {
            radChosen = &shell.radiative.back();
        }

        if (radChosen != nullptr) {
            if (gammaDef != nullptr && radChosen->energy > 0.0) {
                secondaries.emplace_back(gammaDef, position, sampleIsotropicDirection(), radChosen->energy);
                emitted += radChosen->energy;
            }
            if (top < static_cast<int>(vacancies.size())) {
                vacancies[top++] = radChosen->originShell;
            }
            continue;
        }

        // 2. Non-radiative transitions (Auger / Coster-Kronig)
        if (shell.auger.empty()) {
            continue;
        }

        const AtomicTransitionManager::AugerTransition* augChosen = nullptr;
        for (const auto& t : shell.auger) {
            r -= t.probability;
            if (r <= 0.0) {
                augChosen = &t;
                break;
            }
        }
        if (augChosen == nullptr) {
            augChosen = &shell.auger.back();
        }

        if (auger && electronDef != nullptr && augChosen->energy > 0.0) {
            secondaries.emplace_back(electronDef, position, sampleIsotropicDirection(), augChosen->energy);
            emitted += augChosen->energy;
        }
        if (top + 1 < static_cast<int>(vacancies.size())) {
            vacancies[top++] = augChosen->originShell;
            vacancies[top++] = augChosen->augerShell;
        }
    }

    return std::max(0.0, initialBinding - emitted);
}