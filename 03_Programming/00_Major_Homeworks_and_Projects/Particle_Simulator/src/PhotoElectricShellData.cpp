#include "PhotoElectricShellData.h"
#include "AtomicTransitionManager.h"
#include "EndfSection.h"
#include "RandomEngine.h"
#include "Units.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

// ==============================================================================
// Construction (ENDF/B-VIII.0 Photo-Atomic Subshells, MF=23 MT=534..599)
// ==============================================================================
PhotoElectricShellData::PhotoElectricShellData(int z) {
    const std::string path = AtomicTransitionManager::getEndfFileName("photoat", z);
    const std::vector<int> sections = EndfSection::findSections(path, 23);
    tables.reserve(sections.size());

    for (int mt : sections) {
        if (mt < 534 || mt > 599) {
            continue;
        }
        EndfSection section(path, 23, mt);
        static_cast<void>(section.readCont());
        EndfSection::Tab1 tab = section.readTab1();
        if (tab.x.size() < 2 || tab.nbt.empty()) {
            continue;
        }

        ShellTable table;
        table.designator = mt - 533;
        table.bindingEnergy = tab.head.c1 * Units::eV;
        table.energy.reserve(tab.x.size());
        for (double xVal : tab.x) {
            table.energy.push_back(xVal * Units::eV);
        }
        table.sigma = std::move(tab.y);
        table.nbt = std::move(tab.nbt);
        table.interpolation = std::move(tab.interpolation);
        tables.push_back(std::move(table));
    }

    if (tables.empty()) {
        throw std::runtime_error("PhotoElectricShellData: no subshell cross sections in " + path);
    }
}

// ==============================================================================
// Subshell Cross Sections & Sampling
// ==============================================================================
// ENDF interpolation laws: 1 histogram, 2 lin-lin, 3 lin-log, 4 log-lin, 5 log-log
double PhotoElectricShellData::interpolate(long code, double x1, double y1, double x2, double y2, double x) {
    const double dx = x2 - x1;
    if (dx <= 0.0) {
        return y2;
    }

    switch (code) {
        case 1:
            return y1;
        case 3:
            return (x1 > 0.0 && x2 > 0.0) ? y1 + (y2 - y1) * std::log(x / x1) / std::log(x2 / x1) : y1;
        case 4:
            return (y1 > 0.0 && y2 > 0.0) ? y1 * std::pow(y2 / y1, (x - x1) / dx) : y1;
        case 5:
            return (x1 > 0.0 && x2 > 0.0 && y1 > 0.0 && y2 > 0.0)
                   ? y1 * std::pow(y2 / y1, std::log(x / x1) / std::log(x2 / x1))
                   : y1;
        default:
            return y1 + (y2 - y1) * (x - x1) / dx;
    }
}

double PhotoElectricShellData::computeCrossSection(int shell, double e) const {
    const ShellTable& table = tables[shell];
    if (e <= table.bindingEnergy) {
        return 0.0;
    }

    const std::vector<double>& x = table.energy;
    const std::vector<double>& y = table.sigma;
    if (e <= x.front()) {
        return y.front();
    }
    if (e >= x.back()) {
        return y.back();
    }

    const auto upper = std::upper_bound(x.begin(), x.end(), e);
    const std::size_t i = static_cast<std::size_t>(upper - x.begin()) - 1;

    long code = table.interpolation.back();
    if (table.nbt.size() > 1) {
        const long pointIndex = static_cast<long>(i + 2);
        for (std::size_t r = 0; r < table.nbt.size(); ++r) {
            if (pointIndex <= table.nbt[r]) {
                code = table.interpolation[r];
                break;
            }
        }
    }

    return std::max(0.0, interpolate(code, x[i], y[i], x[i + 1], y[i + 1], e));
}

double PhotoElectricShellData::computeShellFraction(int shell, double e) const {
    const int n = getNumberOfShells();
    if (shell < 0 || shell >= n) {
        return 0.0;
    }

    double target = 0.0;
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
        const double s = computeCrossSection(i, e);
        if (i == shell) {
            target = s;
        }
        total += s;
    }
    return total > 0.0 ? target / total : 0.0;
}

int PhotoElectricShellData::sampleShell(double e) const {
    const int n = getNumberOfShells();
    std::array<double, 66> sigma{};
    double total = 0.0;

    for (int i = 0; i < n; ++i) {
        sigma[static_cast<std::size_t>(i)] = computeCrossSection(i, e);
        total += sigma[static_cast<std::size_t>(i)];
    }
    if (total <= 0.0) {
        return -1;
    }

    double r = RandomEngine::uniNumGen() * total;
    for (int i = 0; i < n; ++i) {
        const double s = sigma[static_cast<std::size_t>(i)];
        r -= s;
        if (r <= 0.0 && s > 0.0) {
            return i;
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        if (sigma[static_cast<std::size_t>(i)] > 0.0) {
            return i;
        }
    }
    return -1;
}

// ==============================================================================
// Element Cache
// ==============================================================================
const PhotoElectricShellData& PhotoElectricShellData::get(int z) {
    if (z >= 1 && z <= 100) {
        static std::array<std::unique_ptr<PhotoElectricShellData>, 101> fastData;
        auto& entry = fastData[static_cast<std::size_t>(z)];
        if (!entry) {
            entry.reset(new PhotoElectricShellData(z));
        }
        return *entry;
    }

    static std::map<int, std::unique_ptr<PhotoElectricShellData>> extraData;
    auto& entry = extraData[z];
    if (!entry) {
        entry.reset(new PhotoElectricShellData(z));
    }
    return *entry;
}