// ==============================================================================
// [STAGE 2] Validation: muons in lead against the PDG tables
// ==============================================================================
// PDG "Atomic and Nuclear Properties", muon table for lead (data/muE_lead_Pb.txt,
// https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_lead_Pb.txt).
//
//   Test 1  dE/dx:           computeDEDX vs the table's "Ionization" column
//   Test 2  CSDA range:      our range table vs the table's "CSDA Range"
//   Test 3  Monte Carlo:     full stepping loop, mean track length = CSDA range?
//                            (no fluctuations / Gauss / Landau / Landau + msc)
//   Test 4  Fluctuations:    shape of the loss in a thin layer (peak, width)
//   Test 5  Msc:             angle behind a lead plate vs Highland, step-size independence
//
// Usage:  ./ValidationApp [seed]
// ==============================================================================

#include "Box.h"
#include "Ionization.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "Particle.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Track.h"
#include "Units.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ------------------------------------------------------------------------------
// One row of the PDG table (only the columns we use)
// ------------------------------------------------------------------------------
struct PdgRow {
    double kineticEnergy;  // MeV
    double ionization;     // MeV cm^2/g
    double totalDEDX;      // MeV cm^2/g (ionization + radiative)
    double csdaRange;      // g/cm^2
};

// Numeric rows have 12 numbers. The two special rows ("Minimum ionization",
// "Muon critical energy") have words in them and are skipped.
std::vector<PdgRow> readPdgTable(const std::string& path) {
    std::vector<PdgRow> rows;
    std::ifstream file(path);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream in(line);
        std::vector<double> values;
        std::string token;
        bool numeric = true;
        while (in >> token) {
            try {
                size_t used = 0;
                values.push_back(std::stod(token, &used));
                if (used != token.size()) numeric = false;
            } catch (...) {
                numeric = false;
            }
        }
        if (numeric && values.size() == 12) {
            rows.push_back({values[0], values[2], values[7], values[8]});
        }
    }
    return rows;
}

// Lead, exactly as in the PDG table header (a = 0.0936 as printed there)
std::shared_ptr<Material> makeLead() {
    DensityEffectData d{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    return std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV, 0.5612 * Units::cm, d);
}

double mean(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s / v.size();
}

double rms(const std::vector<double>& v) {
    const double m = mean(v);
    double s = 0.0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / v.size());
}

// Peak (most probable value) and FWHM of a sample, from a histogram
void peakAndFwhm(const std::vector<double>& v, double lo, double hi, int nBins, double& peak, double& fwhm) {
    std::vector<int> h(nBins, 0);
    const double w = (hi - lo) / nBins;
    for (double x : v) {
        const int b = static_cast<int>((x - lo) / w);
        if (b >= 0 && b < nBins) h[b]++;
    }
    // light smoothing (3-bin average) so the peak is not one noisy bin
    std::vector<double> s(nBins, 0.0);
    for (int i = 1; i + 1 < nBins; ++i) s[i] = (h[i - 1] + h[i] + h[i + 1]) / 3.0;
    const int iMax = static_cast<int>(std::max_element(s.begin(), s.end()) - s.begin());
    peak = lo + (iMax + 0.5) * w;
    const double half = 0.5 * s[iMax];
    int iL = iMax, iR = iMax;
    while (iL > 0 && s[iL] > half) --iL;
    while (iR < nBins - 1 && s[iR] > half) ++iR;
    fwhm = (iR - iL) * w;
}

void printHeader(const char* title) {
    std::printf("\n==============================================================================\n");
    std::printf("  %s\n", title);
    std::printf("==============================================================================\n");
}

} // namespace

int main(int argc, char** argv) {
    RandomEngine::setSeed(argc > 1 ? std::stoull(argv[1]) : 12345);
    std::printf("Validation: muons in lead vs PDG.  Seed %llu\n",
                static_cast<unsigned long long>(RandomEngine::getSeed()));

    const std::vector<PdgRow> pdg = readPdgTable(PDG_TABLE_PATH);
    if (pdg.empty()) {
        std::printf("Could not read %s\n", PDG_TABLE_PATH);
        return 1;
    }
    std::printf("Read %zu rows from %s\n", pdg.size(), PDG_TABLE_PATH);

    auto lead = makeLead();
    const double rho = lead->getDensity();
    Particle muon(13, "mu-", PhysicalConstants::m_mu, -1.0, 0.5, -1.0); // stable: Decay is not tested here

    Ionization ionization;
    ionization.setVerboseLevel(0);
    const double unitDEDX = Units::MeV * rho; // GeV/cm  ->  MeV cm^2/g : divide by (MeV * rho)

    // ==========================================================================
    // TEST 1: dE/dx
    // ==========================================================================
    printHeader("TEST 1  dE/dx  [MeV cm^2/g]   ours vs PDG 'Ionization' column");
    std::printf("%12s %10s %10s %9s\n", "T [MeV]", "ours", "PDG", "diff");
    double maxDev = 0.0, maxDevT = 0.0;
    for (const PdgRow& row : pdg) {
        const double ours = ionization.computeDEDX(muon, row.kineticEnergy * Units::MeV, *lead) / unitDEDX;
        const double dev = ours / row.ionization - 1.0;
        const bool dependable = row.kineticEnergy >= 10.0 && row.kineticEnergy <= 1.0e5;
        if (dependable && std::abs(dev) > std::abs(maxDev)) { maxDev = dev; maxDevT = row.kineticEnergy; }

        // print ~3 rows per decade: mantissa 1, 2 or 5
        const double mantissa = row.kineticEnergy / std::pow(10.0, std::floor(std::log10(row.kineticEnergy) + 1e-9));
        const bool show = std::abs(mantissa - 1.0) < 1e-6 || std::abs(mantissa - 2.0) < 1e-6 || std::abs(mantissa - 5.0) < 1e-6;
        if (show && row.kineticEnergy <= 1.0e5) {
            std::printf("%12.4g %10.4f %10.4f %+8.2f%%%s\n", row.kineticEnergy, ours, row.ionization, 100.0 * dev,
                        row.kineticEnergy < 10.0 ? "   (PDG: below 10 MeV not dependable)" : "");
        }
    }
    std::printf("Largest deviation for 10 MeV <= T <= 100 GeV: %+.2f%% at T = %g MeV\n", 100.0 * maxDev, maxDevT);

    // ==========================================================================
    // TEST 2: CSDA range
    // ==========================================================================
    // PDG's CSDA range uses the TOTAL dE/dx (with brems, pair, photonuclear); ours
    // has ionisation only. So we also integrate PDG's ionisation column ourselves
    // ("PDG ion-only"): that is the apples-to-apples comparison. Start value: PDG's
    // range at 10 MeV (radiative losses are 1e-5 of the total there).
    printHeader("TEST 2  CSDA range  [g/cm^2]");
    std::vector<double> pdgIonRange(pdg.size(), 0.0);
    size_t start = 0;
    while (pdg[start].kineticEnergy < 10.0) ++start;
    pdgIonRange[start] = pdg[start].csdaRange;
    for (size_t i = start + 1; i < pdg.size(); ++i) {
        // trapezoid in ln T:  dR = (T/S) dlnT
        const double f0 = pdg[i - 1].kineticEnergy / pdg[i - 1].ionization;
        const double f1 = pdg[i].kineticEnergy / pdg[i].ionization;
        pdgIonRange[i] = pdgIonRange[i - 1] + 0.5 * (f0 + f1) * std::log(pdg[i].kineticEnergy / pdg[i - 1].kineticEnergy);
    }
    std::printf("%12s %11s %11s %9s %13s %9s\n", "T [MeV]", "ours", "PDG ion-only", "diff", "PDG (total)", "diff");
    for (size_t i = start; i < pdg.size(); ++i) {
        const double t = pdg[i].kineticEnergy;
        const double mantissa = t / std::pow(10.0, std::floor(std::log10(t) + 1e-9));
        if (std::abs(mantissa - 1.0) > 1e-6 || t > 1.0e5) continue;
        const double ours = ionization.getRange(muon, t * Units::MeV, *lead) * rho;
        std::printf("%12.4g %11.4g %11.4g %+8.2f%% %13.4g %+8.2f%%\n", t, ours, pdgIonRange[i],
                    100.0 * (ours / pdgIonRange[i] - 1.0), pdg[i].csdaRange, 100.0 * (ours / pdg[i].csdaRange - 1.0));
    }
    std::printf("(The growing difference to 'PDG (total)' above 1 GeV is the missing radiative losses.)\n");

    // ==========================================================================
    // TEST 3: Monte Carlo range with the full stepping loop
    // ==========================================================================
    printHeader("TEST 3  Monte Carlo: muons stopping in a lead block");
    std::printf("Path length L and depth z (projection on the beam axis), in g/cm^2. Expect <L> = CSDA range.\n");
    std::printf("%9s %-14s %7s %11s %11s %8s %9s %11s\n", "T", "physics", "events", "<L>", "CSDA", "diff", "rms(L)/<L>", "<z>/<L>");

    struct Config { const char* name; FluctuationModel model; bool msc; };
    const Config configs[] = {
        {"no fluct", FluctuationModel::none, false},
        {"gauss", FluctuationModel::gauss, false},
        {"landau", FluctuationModel::landau, false},
        {"landau + msc", FluctuationModel::landau, true},
    };
    const double energies[] = {100.0 * Units::MeV, 1.0 * Units::GeV, 10.0 * Units::GeV};

    for (double energy : energies) {
        for (const Config& config : configs) {
            RunManager runManager;
            auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
            runManager.setWorldVolume(std::make_unique<Box>("World", 10.0 * Units::m, 10.0 * Units::m, 10.0 * Units::m, vacuum));
            runManager.addVolume(std::make_unique<Box>("LeadBlock", 9.5 * Units::m, 9.5 * Units::m, 9.5 * Units::m, lead));
            if (config.msc) runManager.addPhysicsProcess(std::make_unique<MultipleScattering>());
            auto ioni = std::make_unique<Ionization>();
            ioni->setFluctuationModel(config.model);
            Ionization* ioniPtr = ioni.get();
            runManager.addPhysicsProcess(std::move(ioni));
            runManager.setVerboseLevel(0);

            const Vector3D startPosition(0.0, 0.0, -9.5 * Units::m); // on the front face: inside (<=)
            double endLength = 0.0, endDepth = 0.0;
            runManager.setUserSteppingAction([&](const Track& track, const Step&) {
                if (!track.isAlive()) {
                    endLength = track.getTrackLength();
                    endDepth = track.getPosition().z - startPosition.z;
                }
            });

            const int nEvents = 4000;
            std::vector<double> lengths, depths;
            for (int i = 0; i < nEvents; ++i) {
                Track primary(&muon, startPosition, Vector3D(0.0, 0.0, 1.0), energy);
                runManager.simulateEvent(primary, 20.0 * Units::m, 100000);
                lengths.push_back(endLength * rho);
                depths.push_back(endDepth * rho);
            }
            const double csda = ioniPtr->getRange(muon, energy, *lead) * rho;
            std::printf("%6g MeV %-14s %7d %11.4g %11.4g %+7.2f%% %9.2f%% %11.4f\n", energy / Units::MeV, config.name,
                        nEvents, mean(lengths), csda, 100.0 * (mean(lengths) / csda - 1.0),
                        100.0 * rms(lengths) / mean(lengths), mean(depths) / mean(lengths));
        }
    }
    std::printf("rms(L)/<L> is the RANGE STRAGGLING (~1-3 %% for muons; 0 without fluctuations).\n");
    std::printf("<z>/<L> < 1 with msc: the zig-zag path is longer than the depth reached.\n");

    // ==========================================================================
    // TEST 4: fluctuation shape in a thin layer
    // ==========================================================================
    printHeader("TEST 4  Energy loss in ONE thin lead layer (sampleFluctuations directly)");
    std::printf("Landau peak vs PDG Delta_p (eq. 34.12) and FWHM vs 4 xi (PDG: w ~ 4 xi). Mean must equal <dE/dx> x.\n");
    std::printf("%8s %7s %-7s %9s %9s %9s %9s %9s %9s %9s\n", "T [GeV]", "x [cm]", "model", "kappa", "mean/<>",
                "peak[MeV]", "Delta_p", "FWHM", "4 xi", "sigma/Bohr");
    const double layerEnergies[] = {1.0 * Units::GeV, 10.0 * Units::GeV};
    const double thicknesses[] = {0.1 * Units::cm, 1.0 * Units::cm};
    for (double energy : layerEnergies) {
        for (double x : thicknesses) {
            const double gamma = 1.0 + energy / muon.getMass();
            const double beta2 = 1.0 - 1.0 / (gamma * gamma);
            const double xi = 0.5 * PhysicalConstants::K * lead->getZoverA() * rho * x / beta2;
            const double tmax = ionization.computeMaxSecondaryEnergy(muon, energy);
            const double meanLoss = ionization.computeDEDX(muon, energy, *lead) * x;
            const double bohrSigma = std::sqrt(xi * tmax * (1.0 - 0.5 * beta2));
            const double deltaP = ionization.computeMostProbableLoss(muon, energy, *lead, x);

            for (FluctuationModel model : {FluctuationModel::gauss, FluctuationModel::landau}) {
                ionization.setFluctuationModel(model);
                std::vector<double> losses;
                const int nSamples = 400000;
                for (int i = 0; i < nSamples; ++i) {
                    losses.push_back(ionization.sampleFluctuations(muon, energy, *lead, meanLoss, x));
                }
                double peak, fwhm;
                peakAndFwhm(losses, 0.0, 3.0 * meanLoss, 600, peak, fwhm);
                std::printf("%8g %7g %-7s %9.2e %9.4f %9.3f %9.3f %9.3f %9.3f %9.3f\n", energy, x,
                            model == FluctuationModel::gauss ? "gauss" : "landau", xi / tmax, mean(losses) / meanLoss,
                            peak / Units::MeV, deltaP / Units::MeV, fwhm / Units::MeV, 4.0 * xi / Units::MeV,
                            rms(losses) / bohrSigma);
            }
        }
    }
    std::printf("For thin layers the Landau peak sits well below the mean, and the Gaussian\n");
    std::printf("(sigma = Bohr's, which is dominated by the rare hard collisions) is far too wide.\n");
    ionization.setFluctuationModel(FluctuationModel::gauss);

    // ==========================================================================
    // TEST 5: multiple scattering behind a lead plate
    // ==========================================================================
    printHeader("TEST 5  Multiple scattering: projected angle behind a lead plate vs Highland");
    std::printf("%8s %8s %8s %12s %12s %12s %8s\n", "T [GeV]", "x [cm]", "x/X0", "max step", "rms(thx)", "Highland", "ratio");
    const double plateEnergies[] = {1.0 * Units::GeV, 10.0 * Units::GeV};
    const double plates[] = {0.1 * Units::cm, 1.0 * Units::cm, 5.0 * Units::cm};
    const double maxSteps[] = {100.0 * Units::cm, 0.01 * Units::cm};
    for (double energy : plateEnergies) {
        for (double x : plates) {
            for (double maxStep : maxSteps) {
                RunManager runManager;
                auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
                runManager.setWorldVolume(std::make_unique<Box>("World", 1.0 * Units::m, 1.0 * Units::m, x / 2.0 + 1.0 * Units::cm, vacuum));
                runManager.addVolume(std::make_unique<Box>("Plate", 0.9 * Units::m, 0.9 * Units::m, x / 2.0, lead));
                runManager.addPhysicsProcess(std::make_unique<MultipleScattering>());
                runManager.addPhysicsProcess(std::make_unique<Ionization>());
                runManager.setVerboseLevel(0);

                double thetaX = 0.0;
                runManager.setUserSteppingAction([&](const Track& track, const Step&) {
                    if (!track.isAlive()) {
                        const Vector3D d = track.getDirection();
                        thetaX = std::atan2(d.x, d.z);
                    }
                });

                const int nEvents = 20000;
                std::vector<double> angles;
                const Vector3D startPosition(0.0, 0.0, -x / 2.0 - 0.5 * Units::cm);
                for (int i = 0; i < nEvents; ++i) {
                    Track primary(&muon, startPosition, Vector3D(0.0, 0.0, 1.0), energy);
                    runManager.simulateEvent(primary, maxStep, 100000);
                    angles.push_back(thetaX);
                }

                Track reference(&muon, startPosition, Vector3D(0.0, 0.0, 1.0), energy);
                MultipleScattering msc;
                const double highland = msc.computeTheta0(reference, *lead, x);
                std::printf("%8g %8g %8.3f %9g cm %9.3f mrad %9.3f mrad %8.3f\n", energy, x, x / lead->getRadLength(),
                            maxStep, rms(angles) * 1e3, highland * 1e3, rms(angles) / highland);
            }
        }
    }
    std::printf("ratio ~1 for BOTH step sizes: the path-length bookkeeping makes the result\n");
    std::printf("independent of the steps. (The muon also loses energy in the plate, so for\n");
    std::printf("thick plates at low energy the ratio is slightly above 1.)\n");

    return 0;
}
