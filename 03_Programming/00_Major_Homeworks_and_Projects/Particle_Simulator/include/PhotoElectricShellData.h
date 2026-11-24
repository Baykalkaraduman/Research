#pragma once

#include <string>
#include <vector>

// ==============================================================================
// Class: PhotoElectricShellData
// Stores ENDF/B-VIII.0 photo-atomic data per element:
//   - MF=23 MT=534+ : subshell photoelectric cross sections [barn vs GeV]
//   - MF=23 MT=502/504 : coherent & incoherent cross sections [barn vs GeV]
//   - MF=27 MT=502/504 : form factor F(x) & scattering function S(x) [vs 1/A]
//   - MF=27 MT=505/506 : imaginary f''(E) & real f'(E) anomalous factors [vs GeV]
// ==============================================================================
class PhotoElectricShellData {
private:
    // --------------------------------------------------------------------------
    // Tabulated Data Structures
    // --------------------------------------------------------------------------
    struct Tab1Table {
        std::vector<double> energy; // [GeV]
        std::vector<double> value;  // [barn] or dimensionless factor
        std::vector<long> nbt;
        std::vector<long> interpolation;
    };

    struct ShellTable {
        int designator = 0;
        double bindingEnergy = 0.0; // [GeV]
        Tab1Table crossSection;
    };

    struct FactorTable {
        std::vector<double> x; // momentum transfer [1/Angstrom]
        std::vector<double> y; // F(x) or S(x)
    };

    std::vector<ShellTable> tables;
    Tab1Table coherent;
    Tab1Table incoherent;
    Tab1Table imaginaryAnomalous; // MF=27, MT=505: f''(E)
    Tab1Table realAnomalous;      // MF=27, MT=506: f'(E)
    FactorTable formFactor;       // MF=27, MT=502: F(x)
    FactorTable scatteringFunction; // MF=27, MT=504: S(x)

    // --------------------------------------------------------------------------
    // Internal Construction & Interpolation
    // --------------------------------------------------------------------------
    explicit PhotoElectricShellData(int z);

    [[nodiscard]] static Tab1Table readEnergyTable(const std::string& path, int mf, int mt,
                                                   bool required, double* c1Out = nullptr);
    [[nodiscard]] static FactorTable readFactorTable(const std::string& path, int mt);
    [[nodiscard]] static double interpolate(long code, double x1, double y1, double x2, double y2, double x);
    [[nodiscard]] static double evaluateTab1(const Tab1Table& table, double e, bool clampNonNegative);
    [[nodiscard]] static double evaluateSubshell(const ShellTable& shell, double e);
    [[nodiscard]] static double interpolateLogLog(const FactorTable& table, double x);

public:
    // --------------------------------------------------------------------------
    // Element Cache & Subshell Accessors
    // --------------------------------------------------------------------------
    [[nodiscard]] static const PhotoElectricShellData& get(int z);

    [[nodiscard]] int getNumberOfShells() const {return static_cast<int>(tables.size());}
    [[nodiscard]] int getDesignator(int shell) const {return tables[static_cast<std::size_t>(shell)].designator;}
    [[nodiscard]] double getBindingEnergy(int shell) const {return tables[static_cast<std::size_t>(shell)].bindingEnergy;}

    // --------------------------------------------------------------------------
    // Shell Probability & Sampling
    // --------------------------------------------------------------------------
    [[nodiscard]] double computeShellFraction(int shell, double e) const;
    [[nodiscard]] int sampleShell(double e) const;

    // --------------------------------------------------------------------------
    // Cross Sections [barn], Form/Scattering Factors, and Anomalous Factors
    // --------------------------------------------------------------------------
    [[nodiscard]] double getCoherentCrossSection(double e) const {return evaluateTab1(coherent, e, true);}
    [[nodiscard]] double getIncoherentCrossSection(double e) const {return evaluateTab1(incoherent, e, true);}
    [[nodiscard]] double getFormFactor(double x) const {return interpolateLogLog(formFactor, x);}
    [[nodiscard]] double getScatteringFunction(double x) const {return interpolateLogLog(scatteringFunction, x);}
    [[nodiscard]] double getImaginaryAnomalousFactor(double e) const {return evaluateTab1(imaginaryAnomalous, e, true);}
    [[nodiscard]] double getRealAnomalousFactor(double e) const {return evaluateTab1(realAnomalous, e, false);}
};