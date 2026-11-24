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
// Construction (ENDF/B-VIII.0 Photo-Atomic: MF=23 Cross Sections & MF=27 Factors)
// ==============================================================================
PhotoElectricShellData::PhotoElectricShellData(int z) {
    const std::string path = AtomicTransitionManager::getEndfFileName("photoat", z);
    const std::vector<int> sections = EndfSection::findSections(path, 23);
    tables.reserve(sections.size());

    for (int mt : sections) {
        if (mt < 534 || mt > 599) {
            continue;
        }
        double bindingEv = 0.0;
        Tab1Table cs = readEnergyTable(path, 23, mt, false, &bindingEv);
        if (!cs.energy.empty()) {
            tables.push_back({mt - 533, bindingEv * Units::eV, std::move(cs)});
        }
    }

    if (tables.empty()) {
        throw std::runtime_error("PhotoElectricShellData: no subshell cross sections in " + path);
    }

    coherent = readEnergyTable(path, 23, 502, true);
    incoherent = readEnergyTable(path, 23, 504, true);
    formFactor = readFactorTable(path, 502);
    scatteringFunction = readFactorTable(path, 504);
    imaginaryAnomalous = readEnergyTable(path, 27, 505, true);
    realAnomalous = readEnergyTable(path, 27, 506, true);
}

PhotoElectricShellData::Tab1Table PhotoElectricShellData::readEnergyTable(const std::string& path, int mf, int mt,
                                                                          bool required, double* c1Out) {
    EndfSection section(path, mf, mt);
    static_cast<void>(section.readCont());
    EndfSection::Tab1 tab = section.readTab1();
    if (tab.x.size() < 2 || tab.x.size() != tab.y.size() || tab.nbt.empty()) {
        if (!required) {
            return {};
        }
        throw std::runtime_error("PhotoElectricShellData: invalid MF=" + std::to_string(mf) +
                                 " MT=" + std::to_string(mt) + " table in " + path);
    }

    if (c1Out != nullptr) {
        *c1Out = tab.head.c1;
    }

    Tab1Table table;
    table.energy.reserve(tab.x.size());
    for (double xVal : tab.x) {
        table.energy.push_back(xVal * Units::eV);
    }
    table.value = std::move(tab.y);
    table.nbt = std::move(tab.nbt);
    table.interpolation = std::move(tab.interpolation);
    return table;
}

PhotoElectricShellData::FactorTable PhotoElectricShellData::readFactorTable(const std::string& path, int mt) {
    EndfSection section(path, 27, mt);
    static_cast<void>(section.readCont());
    EndfSection::Tab1 tab = section.readTab1();
    if (tab.x.size() < 2 || tab.x.size() != tab.y.size()) {
        throw std::runtime_error("PhotoElectricShellData: invalid MF=27 MT=" + std::to_string(mt) + " table in " + path);
    }
    return {std::move(tab.x), std::move(tab.y)};
}

// ==============================================================================
// Interpolation & Table Evaluation
// ==============================================================================
// ENDF interpolation laws: 1 histogram, 2 lin-lin, 3 lin-log, 4 log-lin, 5 log-log
double PhotoElectricShellData::interpolate(long code, double x1, double y1, double x2, double y2, double x) {
    const double dx = x2 - x1;
    if (dx <= 0.0) {
        return y2;
    }

    const double lin = y1 + (y2 - y1) * (x - x1) / dx;
    switch (code) {
        case 1:
            return y1;
        case 3:
            return (x1 > 0.0 && x2 > 0.0) ? y1 + (y2 - y1) * std::log(x / x1) / std::log(x2 / x1) : lin;
        case 4:
            return (y1 > 0.0 && y2 > 0.0) ? y1 * std::pow(y2 / y1, (x - x1) / dx) : lin;
        case 5:
            return (x1 > 0.0 && x2 > 0.0 && y1 > 0.0 && y2 > 0.0)
                   ? y1 * std::pow(y2 / y1, std::log(x / x1) / std::log(x2 / x1))
                   : lin;
        default:
            return lin;
    }
}

double PhotoElectricShellData::evaluateTab1(const Tab1Table& table, double e, bool clampNonNegative) {
    if (table.energy.size() < 2 || (clampNonNegative && e <= 0.0)) {
        return 0.0;
    }

    const std::vector<double>& x = table.energy;
    const std::vector<double>& y = table.value;
    if (e <= x.front()) {
        return clampNonNegative ? std::max(0.0, y.front()) : y.front();
    }
    if (e >= x.back()) {
        return clampNonNegative ? std::max(0.0, y.back()) : y.back();
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

    const double val = interpolate(code, x[i], y[i], x[i + 1], y[i + 1], e);
    return clampNonNegative ? std::max(0.0, val) : val;
}

double PhotoElectricShellData::evaluateSubshell(const ShellTable& shell, double e) {
    if (e <= shell.bindingEnergy) {
        return 0.0;
    }
    return evaluateTab1(shell.crossSection, e, true);
}

double PhotoElectricShellData::interpolateLogLog(const FactorTable& table, double x) {
    if (table.x.empty()) {
        return 0.0;
    }
    if (x <= table.x.front()) {
        return table.y.front();
    }
    if (x >= table.x.back()) {
        return table.y.back();
    }

    const auto upper = std::upper_bound(table.x.begin(), table.x.end(), x);
    const std::size_t i = static_cast<std::size_t>(upper - table.x.begin()) - 1;
    const double x1 = table.x[i];
    const double x2 = table.x[i + 1];
    const double y1 = table.y[i];
    const double y2 = table.y[i + 1];
    const double dx = x2 - x1;
    if (dx <= 0.0) {
        return y2;
    }

    if (x1 <= 0.0 || x2 <= 0.0 || y1 <= 0.0 || y2 <= 0.0) {
        return y1 + (y2 - y1) * (x - x1) / dx;
    }

    const double slope = std::log(y2 / y1) / std::log(x2 / x1);
    return y1 * std::pow(x / x1, slope);
}

// ==============================================================================
// Subshell Fractions & Sampling
// ==============================================================================
double PhotoElectricShellData::computeShellFraction(int shell, double e) const {
    const int n = getNumberOfShells();
    if (shell < 0 || shell >= n) {
        return 0.0;
    }

    double target = 0.0;
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
        const double s = evaluateSubshell(tables[static_cast<std::size_t>(i)], e);
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
        const double s = evaluateSubshell(tables[static_cast<std::size_t>(i)], e);
        sigma[static_cast<std::size_t>(i)] = s;
        total += s;
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