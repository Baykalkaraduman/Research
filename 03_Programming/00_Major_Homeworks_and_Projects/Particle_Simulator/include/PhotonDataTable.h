#pragma once

#include <vector>

// ==============================================================================
// Class: PhotonDataTable
// Tabulated photon cross sections of one element, read from
// data/xcom_Z<Z>.txt  (NIST XCOM: E, coherent, incoherent, photoelectric, pair, total)
//
// Used for:
//   photoelectric: no formula covers the absorption edges (PhotoElectricEffect)
//   pair production: the Bethe-Heitler formula (Born approximation) is 30-80 %
//     too low near threshold and ~10 % too high at 10-50 MeV in lead, because
//     the nuclear field is strong (alpha Z = 0.6). Geant4's G4BetheHeitlerModel
//     takes its cross section from an empirical fit to such data, too.
// The Compton cross section is NOT taken from here: it is Klein-Nishina.
//
// One table per Z, read once and kept for the whole run (get()).
// ==============================================================================
class PhotonDataTable {
private:
    std::vector<double> energy;        // GeV, ascending; an absorption edge appears twice
    std::vector<double> photoelectric; // mu/rho, cm^2/g
    std::vector<double> pair;          // mu/rho, cm^2/g (nuclear + electron field)
    std::vector<double> shellEnergies; // the edge energies = binding energies, ascending

    explicit PhotonDataTable(int z); // reads the file (-> PhotonDataTable.cpp)

    // log-log interpolation in one column (-> PhotonDataTable.cpp)
    double interpolate(const std::vector<double>& values, double e) const;

public:
    // The table of element Z, read on first use
    static const PhotonDataTable& get(int z);

    // mu/rho in cm^2/g
    double getPhotoElectric(double e) const { return interpolate(photoelectric, e); }
    double getPair(double e) const { return interpolate(pair, e); }

    const std::vector<double>& getShellEnergies() const { return shellEnergies; }
};
