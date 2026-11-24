#include "PhotonDataTable.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#ifndef TOYSIM_DATA_DIR
#define TOYSIM_DATA_DIR "data"
#endif

// ==============================================================================
// Internal Helpers
// ==============================================================================
namespace {

// Parses the first 5 numeric columns of an XCOM row without allocating a stream.
bool parseXcomRow(std::string_view line, double& energyMeV, double& pe, double& pp) noexcept {
    const char* ptr = line.data();
    const char* end = ptr + line.size();

    while (ptr < end && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r')) {
        ++ptr;
    }
    if (ptr == end || *ptr == '#') {
        return false;
    }

    double cols[5] = {};
    for (int c = 0; c < 5; ++c) {
        char* next = nullptr;
        cols[c] = std::strtod(ptr, &next);
        if (next == ptr || next > end) {
            return false;
        }
        ptr = next;
    }

    energyMeV = cols[0];
    pe = cols[3];
    pp = cols[4];
    return true;
}

} // namespace

// ==============================================================================
// Construction (NIST XCOM Cross-Section Table)
// ==============================================================================
PhotonDataTable::PhotonDataTable(int z) {
    char path[256];
    std::snprintf(path, sizeof(path), "%s/xcom_Z%d.txt", TOYSIM_DATA_DIR, z);

    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(std::string("PhotonDataTable: no cross-section table ") + path);
    }

    energy.reserve(128);
    photoelectric.reserve(128);
    pair.reserve(128);

    std::string line;
    double energyMeV = 0.0;
    double pe = 0.0;
    double pp = 0.0;

    while (std::getline(in, line)) {
        if (!parseXcomRow(line, energyMeV, pe, pp)) {
            continue;
        }
        const double e = energyMeV * Units::MeV;
        if (!energy.empty() && e == energy.back()) {
            shellEnergies.push_back(e);
        }
        energy.push_back(e);
        photoelectric.push_back(pe);
        pair.push_back(pp);
    }

    if (energy.size() < 2) {
        throw std::runtime_error(std::string("PhotonDataTable: table ") + path + " is empty");
    }
}

// ==============================================================================
// Element Cache
// ==============================================================================
const PhotonDataTable& PhotonDataTable::get(int z) {
    if (z >= 1 && z <= 100) {
        static std::array<std::unique_ptr<PhotonDataTable>, 101> fastTables;
        auto& entry = fastTables[static_cast<std::size_t>(z)];
        if (!entry) {
            entry.reset(new PhotonDataTable(z));
        }
        return *entry;
    }

    static std::map<int, std::unique_ptr<PhotonDataTable>> extraTables;
    auto& entry = extraTables[z];
    if (!entry) {
        entry.reset(new PhotonDataTable(z));
    }
    return *entry;
}

// ==============================================================================
// Log-Log Power-Law Interpolation
// ==============================================================================
double PhotonDataTable::interpolate(const std::vector<double>& values, double e) const {
    if (e <= energy.front()) {
        return values.front();
    }

    std::size_t j = static_cast<std::size_t>(std::upper_bound(energy.begin(), energy.end(), e) - energy.begin());
    if (j >= energy.size()) {
        j = energy.size() - 1;
    }
    const std::size_t i = j - 1;
    const double de = energy[j] - energy[i];
    if (de <= 0.0) {
        return values[j];
    }

    if (values[i] <= 0.0 || values[j] <= 0.0) {
        return values[i] + (values[j] - values[i]) * (e - energy[i]) / de;
    }

    const double slope = std::log(values[j] / values[i]) / std::log(energy[j] / energy[i]);
    return values[i] * std::pow(e / energy[i], slope);
}