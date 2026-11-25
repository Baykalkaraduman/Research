#pragma once

#include "Material.h"

// ==============================================================================
// Utility: LPMFunctions
// Migdal (1956) & Klein (1999) suppression factors for high-energy
// bremsstrahlung (LPM and/or Ter-Mikaelian dielectric effect) and pair production.
// Geant4's xi*phi <= 1 clamp guarantees R <= 1 everywhere (no rejection envelope needed).
// ==============================================================================
class LPMFunctions {
public:
    // Characteristic scales [GeV]
    [[nodiscard]] static double computeLPMEnergy(const Material& material);
    [[nodiscard]] static double computePlasmaEnergy(const Material& material);

    // Migdal functions: xi(s) in [1, 2], phi(s) in [0, 1], G(s) in [0, 1]
    [[nodiscard]] static double computeS1(int z);
    [[nodiscard]] static double computeXi(double s, double s1);
    [[nodiscard]] static double solveMigdalS(double sPrime, double s1);
    [[nodiscard]] static double computePhi(double s);
    [[nodiscard]] static double computeG(double s);

    // Suppression ratios relative to Bethe-Heitler (guaranteed in [0, 1])
    [[nodiscard]] static double computeDielectricFactor(double electronTotalEnergy, double photonEnergy,
                                                        const Material& material);
    [[nodiscard]] static double computeBremsstrahlungFactor(double electronTotalEnergy, double photonEnergy,
                                                            const Material& material, bool includeDielectric = true);
    [[nodiscard]] static double computePairProductionFactor(double photonEnergy, double pairFraction,
                                                            const Material& material);
};