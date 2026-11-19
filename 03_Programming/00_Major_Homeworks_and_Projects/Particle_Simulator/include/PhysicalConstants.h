#pragma once

#include "Units.h"
namespace PhysicalConstants {

    // Speed of Light
    constexpr double c = 29.9792458 * Units::cm / Units::ns;

    // Mass of Electron
    constexpr double m_e = 0.51099895 * Units::MeV;

    // ==================================
    // Muon Lepton
    // ==================================
    // Mass of Muon
    constexpr double m_mu = 105.6583755 * Units::MeV;
    // Mean Lifetime of Muon
    constexpr double tau_mu = 2.1969811 * Units::us;

    // ==================================
    // Tau Lepton
    // ==================================
    // Mass of Tau
    constexpr double m_tau = 1.77686 * Units::GeV;
    // Mean Lifetime of Tau
    constexpr double tau_tau = 290.3 * Units::fs;

    // Bethe-Bloch prefactor K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm^2/mol
    // density is always in g/cm^3 and A in g/mol. Then K * (Z/A) * density has
    // the unit MeV/cm, exactly what dE/dx needs.
    constexpr double K = 0.307075 * Units::MeV;

    // Constants of electromagnetic interactions (CODATA 2018)
    // Every photon / brems / pair cross section is a number times r_e^2
    // (and alpha for the radiative ones).
    constexpr double alpha = 1.0 / 137.035999084;                   // fine-structure constant
    constexpr double r_e = 2.8179403262e-13 * Units::cm;            // classical electron radius
    constexpr double N_A = 6.02214076e23;                           // Avogadro (1/mol)

    // Highland formula constant (PDG eq. 34.15)
    constexpr double highlandConstant = 13.6 * Units::MeV;
}