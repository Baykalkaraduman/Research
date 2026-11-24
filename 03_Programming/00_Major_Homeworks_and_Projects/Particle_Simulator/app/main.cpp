#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "AnalysisManager.h"
#include "Annihilation.h"
#include "AtomicDeexcitation.h"
#include "Box.h"
#include "Bremsstrahlung.h"
#include "ComptonScattering.h"
#include "Decay.h"
#include "DecayTable.h"
#include "GammaConversion.h"
#include "Ionization.h"
#include "Material.h"
#include "MultipleScattering.h"
#include "Particle.h"
#include "ParticleTable.h"
#include "PhotoElectricEffect.h"
#include "PhysicalConstants.h"
#include "RandomEngine.h"
#include "RunManager.h"
#include "Track.h"
#include "Units.h"
#include "Vector3D.h"

int main(int argc, char** argv) {
    std::cout << "========================================================\n"
              << "   TOY Particle Simulation\n"
              << "========================================================\n\n";

    // ==============================================================================
    // 0. Random Seed
    // ==============================================================================
    const RandomEngine::seedType seed = (argc > 1) ? std::stoull(argv[1]) : std::random_device{}();
    RandomEngine::setSeed(seed);
    std::cout << "[RandomEngine] Seed: " << RandomEngine::getSeed() << "\n\n";

    // ==============================================================================
    // 1. Materials & Geometry
    // ==============================================================================
    // Lead (PDG parameters; X0 = 6.37 g/cm^2 / 11.35 g/cm^3 = 0.5612 cm)
    const DensityEffectData leadDensityEffect{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    auto lead = std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV,
                                           0.5612 * Units::cm, leadDensityEffect);
    auto vacuum = std::make_shared<Material>("Vacuum", 0.0);

    RunManager runManager;
    runManager.setWorldVolume(std::make_unique<Box>("World", 1.0 * Units::m, 1.0 * Units::m, 1.0 * Units::m, vacuum));
    runManager.addVolume(std::make_unique<Box>("Target_LeadBlock", 50.0, 50.0, 10.0, lead)); // Z in [-10, +10] cm

    // ==============================================================================
    // 2. Particle Definitions
    // ==============================================================================
    ParticleTable particleTable;
    const Particle& electronDef = *particleTable.findParticle("e-");
    const Particle& positronDef = *particleTable.findParticle("e+");
    const Particle& gammaDef = *particleTable.findParticle("gamma");

    // Short-lived Toy Muon (tau = 22 ps, c*tau = 6.6 mm) using real mu- mass and Michel decay table
    Particle toyMuonDef(13, "Toy_Muon", PhysicalConstants::m_mu, -1.0, 0.5, 0.022);
    toyMuonDef.setDecayTable(particleTable.findParticle("mu-")->getDecayTable());

    // ==============================================================================
    // 3. Physics Processes
    // ==============================================================================
    // MultipleScattering must precede Ionization (uses pre-step momentum)
    runManager.addPhysicsProcess(std::make_unique<MultipleScattering>());

    auto ionization = std::make_unique<Ionization>();
    ionization->setFluctuationModel(FluctuationModel::landau);
    ionization->setDeltaRayProduction(&electronDef, 1.0 * Units::MeV);
    runManager.addPhysicsProcess(std::move(ionization));

    runManager.addPhysicsProcess(std::make_unique<Decay>());
    runManager.addPhysicsProcess(std::make_unique<Bremsstrahlung>(&gammaDef));
    runManager.addPhysicsProcess(std::make_unique<Annihilation>(&gammaDef));

    // Photoelectric absorption with atomic relaxation (fluorescence on, Auger off by default)
    AtomicDeexcitation deexcitation(&gammaDef, &electronDef);
    auto photoElectric = std::make_unique<PhotoElectricEffect>(&electronDef);
    photoElectric->setAtomicDeexcitation(&deexcitation);
    runManager.addPhysicsProcess(std::move(photoElectric));

    runManager.addPhysicsProcess(std::make_unique<ComptonScattering>(&electronDef));
    runManager.addPhysicsProcess(std::make_unique<GammaConversion>(&electronDef, &positronDef));

    // ==============================================================================
    // 4. Primary Track Initialization
    // ==============================================================================
    const std::string primaryName = (argc > 2) ? argv[2] : "mu-";
    const double primaryEnergy = (argc > 3) ? std::stod(argv[3]) * Units::GeV : 5.0 * Units::GeV;
    const int verbose = (argc > 4) ? std::stoi(argv[4]) : 1;

    const Particle* primaryDef = &toyMuonDef;
    if (primaryName != "mu-" && primaryName != "Toy_Muon") {
        primaryDef = particleTable.findParticle(primaryName);
        if (primaryDef == nullptr) {
            std::cerr << "Unknown particle: " << primaryName << "\n";
            return 1;
        }
    }
    runManager.setVerboseLevel(verbose);

    std::cout << ">>> INJECTING PRIMARY " << primaryDef->getName()
              << " (" << primaryEnergy << " GeV) INTO THE VACUUM >>>\n";

    const Vector3D initialPosition(0.0, 0.0, -20.0);
    const Vector3D initialDirection(0.0, 0.0, 1.0);
    Track primaryTrack(primaryDef, initialPosition, initialDirection, primaryEnergy);

    // ==============================================================================
    // 5. Data Analysis Pipeline
    // ==============================================================================
    AnalysisManager analysis("simulation_output.root");
    analysis.initializeHistograms();
    runManager.setAnalysisManager(&analysis);

    // ==============================================================================
    // 6. Event Execution & Energy Bookkeeping
    // ==============================================================================
    const double maxStepLength = 10.0 * Units::cm;
    const int maxStepsPerTrack = 1000;

    // Balance: deposited + escaped + neutrinos = T0 + Q(decays) + 2*m_e*(N_annih - N_conv)
    double deposited = 0.0;
    double escaped = 0.0;
    double escapedNeutrino = 0.0;
    double decayEnergy = 0.0;

    runManager.setUserSteppingAction([&](const Track& track, const Step& step) {
        deposited += step.totalEnergyDeposit;
        if (step.leftWorld) {
            const int pdg = std::abs(track.getParticleDefinition()->getPdgId());
            if (pdg == 12 || pdg == 14 || pdg == 16) {
                escapedNeutrino += track.getKineticEnergy();
            } else {
                escaped += track.getKineticEnergy();
            }
        }
        if (step.processDefinedStep == "Decay" && step.secondaryInCurrentStep != nullptr) {
            double q = track.getMass();
            for (const Track& d : *step.secondaryInCurrentStep) {
                q -= d.getMass();
            }
            decayEnergy += q;
        }
    });

    runManager.simulateEvent(primaryTrack, maxStepLength, maxStepsPerTrack);

    std::cout << "\n[Summary] tracks: " << runManager.getNumberOfTracks()
              << " | largest stack: " << runManager.getMaxStackSize()
              << " | deposited: " << deposited << " GeV | escaped: " << escaped << " GeV"
              << " | neutrinos: " << escapedNeutrino << " GeV | released by decays: " << decayEnergy << " GeV\n";

    // ==============================================================================
    // 7. Save Output
    // ==============================================================================
    analysis.saveAndClose();
    std::cout << "Simulation finished successfully.\n";
    return 0;
}