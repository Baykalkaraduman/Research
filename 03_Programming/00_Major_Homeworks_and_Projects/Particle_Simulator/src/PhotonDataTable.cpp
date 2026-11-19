#include "PhotonDataTable.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

// Folder of the data tables, set by CMake
#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

PhotonDataTable::PhotonDataTable(int z) {
    const std::string path = std::string(TOYSIM_DATA_DIR) + "/xcom_Z" + std::to_string(z) + ".txt";
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Photon Data Table: no cross-section table " + path);
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream columns(line);
        double energyMeV, coherent, incoherent, pe, pp;
        if (!(columns >> energyMeV >> coherent >> incoherent >> pe >> pp)) continue;
        const double e = energyMeV * Units::MeV;
        if (!energy.empty() && e == energy.back()) {
            shellEnergies.push_back(e); // second row of an absorption edge
        }
        energy.push_back(e);
        photoelectric.push_back(pe);
        pair.push_back(pp);
    }
    if (energy.size() < 2) {
        throw std::runtime_error("Photon Data Table: table " + path + " is empty");
    }
}

// One table per Z for the whole run. The function-local static map is created
// at the first call; each element's file is read the first time it is asked for.
const PhotonDataTable& PhotonDataTable::get(int z) {
    static std::map<int, std::unique_ptr<PhotonDataTable>> tables;
    auto found = tables.find(z);
    if (found == tables.end()) {
        found = tables.emplace(z, std::unique_ptr<PhotonDataTable>(new PhotonDataTable(z))).first;
    }
    return *found->second;
}

// Cross sections follow power laws, so straight lines in log(mu) vs log(E) (what
// XCOM itself uses). upper_bound returns the first row with a HIGHER energy; at
// an edge both rows are skipped, so an energy exactly on the edge already sees
// the value above it. Where a value is 0 (pair production up to 1.022 MeV)
// there is no logarithm: linear interpolation there.
double PhotonDataTable::interpolate(const std::vector<double>& values, double e) const {
    if (e <= energy.front()) {
        return values.front();
    }
    std::size_t j = std::upper_bound(energy.begin(), energy.end(), e) - energy.begin();
    if (j >= energy.size()) {
        j = energy.size() - 1; // above the table: extend the last power law
    }
    const std::size_t i = j - 1;
    if (values[i] <= 0.0 || values[j] <= 0.0) {
        return values[i] + (values[j] - values[i]) * (e - energy[i]) / (energy[j] - energy[i]);
    }
    const double slope = std::log(values[j] / values[i]) / std::log(energy[j] / energy[i]);
    return values[i] * std::pow(e / energy[i], slope);
}
