#pragma once

#include "GeoDefs.h"
#include "PhysicalConstants.h"

#include <string>
#include <iostream>

// ==============================================================================
// Class: Material
// Represents the physical medium a particle travels through.
// Equivalent to G4Material.
// ==============================================================================
struct DensityEffectData {
    double cBar   = 0.0;   // -C  (the table prints it as "Cbar")
    double x0     = 0.0;   // below x0: (almost) no density effect
    double x1     = 0.0;   // above x1: the high-energy asymptote 2 ln10 x - cBar
    double a      = 0.0;
    double k      = 0.0;   // (the table prints it as "k=m_s")
    double delta0 = 0.0;   // conductors only: delta at x0 (0 for insulators)
};

class Material {
private:
    std::string name;
    double density;               // g/cm^3
    double Z;                     // atomic number Z
    double A;                     // molar mass A in g/mol
    double meanExcitationEnergy;  // I (G4IonisParamMat::GetMeanExcitationEnergy)
    double radLength;             // radiation length X0 in cm (G4Material::GetRadlen)
    DensityEffectData densityEffect;

public:
    // All-arg Constructor
    Material(std::string matName, double matDensity, double matZ = 0.0, double matA = 1.0,
            double matI = 0.0, double matRadLength = GeoDefs::kInfinity, DensityEffectData matDensityEffect = {})
            : name(std::move(matName)), density(matDensity), Z(matZ), A(matA),
            meanExcitationEnergy(matI), radLength(matRadLength), densityEffect(matDensityEffect) {}
    // ========================================
    // Getters
    // ========================================

    std::string getName() const     {return name;}
    double getDensity() const       {return density;}
    double getZ() const             {return Z;}
    double getA() const             {return A;}
    double getMeanExcEn() const&    {return meanExcitationEnergy;}
    double getRadLength() const&    {return radLength;}
    const DensityEffectData getDensEffData() const {return densityEffect;}

    // (Z / A) [mol/g] to use in Bethe-Bloch formula
    double getZoverA() const {return Z/A;}

    // A cross section per atom (cm^2) times n is the inverse mean free path.
    // Atoms per cm^3: n = rho N_A / A    
    double getAtomDensity() const { return A > 0.0 ? density * PhysicalConstants::N_A / A : 0.0; }

    // The Born approximation behind the Bethe-Heitler formulas treats the
    // nuclear field as weak; for lead (alpha Z = 0.6) it is not, and pair
    // production / brems at high energy are ~9 % smaller.
    double computeCoulombCorrection() const;

    // ========================================
    // Functions
    // ========================================
    // Density-effect correction delta at x = log10(beta*gamma)
    double DensityCorrection(double x) const;
    
    // ========================================
    // Utility
    // ========================================
    void display() const {
        std::cout << "Material: " << name << " [Density: " << density << " g/cm^3]\n";
    }
};