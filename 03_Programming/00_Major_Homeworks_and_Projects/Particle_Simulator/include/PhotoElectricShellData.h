#pragma once

#include <vector>

// ==============================================================================
// Class: PhotoElectricShellData
// Stores ENDF/B-VIII.0 (MF=23, MT=534+) subshell photoelectric cross sections
// per element to sample which atomic shell absorbs a photon.
// ==============================================================================
class PhotoElectricShellData {
private:
    // --------------------------------------------------------------------------
    // Subshell Table & Member Data
    // --------------------------------------------------------------------------
    struct ShellTable {
        int designator = 0;
        double bindingEnergy = 0.0; // [GeV]
        std::vector<double> energy; // [GeV]
        std::vector<double> sigma;  // [barn]
        std::vector<long> nbt;
        std::vector<long> interpolation;
    };

    std::vector<ShellTable> tables;

    // --------------------------------------------------------------------------
    // Internal Construction & Interpolation
    // --------------------------------------------------------------------------
    explicit PhotoElectricShellData(int z);

    [[nodiscard]] static double interpolate(long code, double x1, double y1, double x2, double y2, double x);
    [[nodiscard]] double computeCrossSection(int shell, double e) const;

public:
    // --------------------------------------------------------------------------
    // Element Cache & Subshell Accessors
    // --------------------------------------------------------------------------
    [[nodiscard]] static const PhotoElectricShellData& get(int z);

    [[nodiscard]] int getNumberOfShells() const {return static_cast<int>(tables.size());}
    [[nodiscard]] int getDesignator(int shell) const {return tables[shell].designator;}
    [[nodiscard]] double getBindingEnergy(int shell) const {return tables[shell].bindingEnergy;}

    // --------------------------------------------------------------------------
    // Shell Probability & Sampling
    // --------------------------------------------------------------------------
    [[nodiscard]] double computeShellFraction(int shell, double e) const;
    [[nodiscard]] int sampleShell(double e) const;
};