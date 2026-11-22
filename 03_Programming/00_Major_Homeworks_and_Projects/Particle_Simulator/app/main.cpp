#include <iostream>
#include <memory>

// Include our custom engine headers
#include "Vector3D.h"
#include "Particle.h"
#include "Track.h"
#include "Material.h"
#include "Box.h"
#include "Ionization.h"
#include "Decay.h"
#include "ParticleTable.h"
#include "DecayTable.h"
#include "MultipleScattering.h"
#include "Bremsstrahlung.h"
#include "Annihilation.h"
#include "PhotoElectricEffect.h"
#include "ComptonScattering.h"
#include "GammaConversion.h"
#include "PhysicalConstants.h"
#include "RunManager.h"
#include "AnalysisManager.h"
#include "Units.h"
#include "RandomEngine.h"

#include <cstdlib>
#include <random>
#include <string>

int main(int argc, char** argv) {
    std::cout << "========================================================\n";
    std::cout << "   TOY Particle Simulation      \n";
    std::cout << "========================================================\n\n";

    // ==========================================================================
    // 0. RANDOM SEED (set ONCE, before anything draws random numbers)
    // ==========================================================================
    RandomEngine::seedType seed;
    if (argc > 1) {
        seed = std::stoull(argv[1]);
    } else {
        seed = std::random_device{}();
    }
    RandomEngine::setSeed(seed);
    std::cout << "[RandomEngine] Seed: " << RandomEngine::getSeed() << "\n\n";

    // ==========================================================================
    // 1. INITIALIZE MATERIALS & GEOMETRY
    // ==========================================================================
    // Lead with the atomic data Bethe-Bloch and Highland need
    // (PDG: header of data/muE_lead_Pb.txt; X0 = 6.37 g/cm^2 / 11.35 g/cm^3)
    DensityEffectData leadDensityEffect{6.2018, 0.3776, 3.8073, 0.0936, 3.1608, 0.14};
    auto lead = std::make_shared<Material>("Lead (Pb)", 11.35, 82.0, 207.2, 823.0 * Units::eV,
                                           0.5612 * Units::cm, leadDensityEffect);
    auto vacuum = std::make_shared<Material>("Vacuum", 0.0);
    RunManager runManager;

    // Worlds volume, Everything happens inside it; a particle that leaves it is killed. 2 m x 2 m x 2 m, filled with vacuum.
    auto world = std::make_unique<Box>("World", 1.0 * Units::m, 1.0 * Units::m, 1.0 * Units::m, vacuum);
    runManager.setWorldVolume(std::move(world));

    // Create a Lead Box. 
    // hZ = 10.0 means it spans from Z = -10.0 to Z = +10.0
    auto leadBlock = std::make_unique<Box>("Target_LeadBlock", 50.0, 50.0, 10.0, lead);
    runManager.addVolume(std::move(leadBlock));

    // ==========================================================================
    // 2. DEFINE STATIC PARTICLE PROPERTIES
    // ==========================================================================
    // All particles and their decay tables come from the ParticleTable
    // (e, gamma, mu, tau, pi, the six neutrinos). Before: each one was
    // defined here by hand, and the muon "decayed" into one electron.
    ParticleTable particleTable;
    const Particle& electronDef = *particleTable.findParticle("e-");
    const Particle& positronDef = *particleTable.findParticle("e+");
    const Particle& gammaDef = *particleTable.findParticle("gamma");

    // We use the short-lived Toy Muon again so you can guarantee seeing the
    // decay happen shortly after it hits the lead block.
    // Real muon mass (was 0.105658), and the real mu- decay table: it now
    // decays into e- anti_nu_e nu_mu with the Michel spectrum. Only the
    // lifetime is "toy": 22 ps instead of 2.2 us (c tau = 6.6 mm, not 659 m).
    Particle toyMuonDef(13, "Toy_Muon", PhysicalConstants::m_mu, -1.0, 0.5, 0.022);
    toyMuonDef.setDecayTable(particleTable.findParticle("mu-")->getDecayTable());

    // ==========================================================================
    // 3. REGISTER PHYSICS PROCESSES
    // ==========================================================================
    // msc BEFORE Ionization: it needs the pre-step momentum (see MultipleScattering.h)
    runManager.addPhysicsProcess(std::make_unique<MultipleScattering>());
    auto ionization = std::make_unique<Ionization>();
    ionization->setFluctuationModel(FluctuationModel::landau);
    // Delta rays above 1 MeV become real electrons (what Geant4's default 0.7 mm range cut gives for e- in lead). Below 1 MeV: continuous loss.
    ionization->setDeltaRayProduction(&electronDef, 1.0 * Units::MeV);
    runManager.addPhysicsProcess(std::move(ionization));

    runManager.addPhysicsProcess(std::move(ionization));
    // No product argument: each particle's DecayTable says what it decays into
    runManager.addPhysicsProcess(std::make_unique<Decay>());

    // Electromagnetic physics of e-, e+ and gamma. Each process says itself which particles it acts on (isApplicable).
    runManager.addPhysicsProcess(std::make_unique<Bremsstrahlung>(&gammaDef));      // e-, e+
    runManager.addPhysicsProcess(std::make_unique<Annihilation>(&gammaDef));        // e+
    runManager.addPhysicsProcess(std::make_unique<PhotoElectricEffect>(&electronDef)); // gamma
    runManager.addPhysicsProcess(std::make_unique<ComptonScattering>(&electronDef));   // gamma
    runManager.addPhysicsProcess(std::make_unique<GammaConversion>(&electronDef, &positronDef)); // gamma

    // e-, e+ and gamma below 100 keV are not followed: their energy stays where they are (a positron still annihilates)
    runManager.setEnergyCut(100.0 * Units::keV);

    // Primary from the command line
    const std::string primaryName = argc > 2 ? argv[2] : "mu-";
    const double primaryEnergy = argc > 3 ? std::stod(argv[3]) * Units::GeV : 5.0 * Units::GeV;
    const int verbose = argc > 4 ? std::stoi(argv[4]) : 1;
    // Any particle of the table; "mu-" stays the toy muon (as before)
    const Particle* primaryDef = &toyMuonDef;
    if (primaryName != "mu-" && primaryName != "Toy_Muon") {
        primaryDef = particleTable.findParticle(primaryName);
        if (primaryDef == nullptr) {
            std::cerr << "Unknown particle: " << primaryName << "\n";
            return 1;
        }
    }
    runManager.setVerboseLevel(verbose);

    // ==========================================================================
    // 4. INITIALIZE THE PRIMARY TRACK
    // ==========================================================================
    std::cout << ">>> INJECTING PRIMARY " << primaryDef->getName() << " (" << primaryEnergy << " GeV) INTO THE VACUUM >>>\n";
    
    // Start particle at Z = -20 (10 cm OUTSIDE the lead block)
    Vector3D initialPosition(0.0, 0.0, -20.0);       
    Vector3D initialDirection(0.0, 0.0, 1.0); // renamed: it is only a direction
    Track primaryTrack(primaryDef, initialPosition, initialDirection, primaryEnergy);

    // =================================================================
    // 5. DATA ANALYSIS PIPELINE SETUP
    // =================================================================
    // Create the Analysis Manager and set the output filename
    AnalysisManager analysis("simulation_output.root");
    
    // Initialize the ROOT histograms in memory
    analysis.initializeHistograms();

    // Inject the Analysis Manager into the Run Manager BEFORE the run starts
    runManager.setAnalysisManager(&analysis);

    // ==========================================================================
    // 6. EXECUTE THE SIMULATION
    // ==========================================================================
    double maxStepLength = 10.0 * Units::cm;

    int maxStepsPerTrack = 1000;

    // Now the RunManager knows about the AnalysisManager and will record data!
    // Energy bookkeeping of the event with a stepping action  + the energy carried away by neutrinos, and the energy the decays
    // release (Q = rest mass turned into kinetic energy). Balance: deposited + escaped + neutrinos = T0 + Q(decays) + 2 m_e (N_annih - N_conv)
    double deposited = 0.0, escaped = 0.0, escapedNeutrino = 0.0, decayEnergy = 0.0;
    runManager.setUserSteppingAction([&](const Track& track, const Step& step) {
        deposited += step.totalEnergyDeposit;
        if (step.leftWorld) {
            const int pdg = std::abs(track.getParticleDefinition()->getPdgId());
            if (pdg == 12 || pdg == 14 || pdg == 16) escapedNeutrino += track.getKineticEnergy();
            else escaped += track.getKineticEnergy();
        }
        if (step.processDefinedStep == "Decay" && step.secondaryInCurrentStep != nullptr) {
            double q = track.getMass();
            for (const Track& d : *step.secondaryInCurrentStep) q -= d.getMass();
            decayEnergy += q;
        }
    });

    // A shower has thousands of tracks: 1000 steps per TRACK is still enough
    runManager.simulateEvent(primaryTrack, maxStepLength, maxStepsPerTrack);

    std::cout << "\n[Summary] tracks: " << runManager.getNumberOfTracks()
              << " | largest stack: " << runManager.getMaxStackSize()
              << " | deposited: " << deposited << " GeV | escaped: " << escaped << " GeV"
              << " | neutrinos: " << escapedNeutrino << " GeV | released by decays: " << decayEnergy << " GeV\n";

    // =================================================================
    // 7. SAVE DATA
    // =================================================================
    // CRITICAL: Write the RAM buffer to the actual .root file on disk
    analysis.saveAndClose();

    std::cout << "Simulation finished successfully.\n";
    return 0;
}