# ==============================================================================
# ToyMonteCarlo — CMake Build Configuration
# ==============================================================================
cmake_minimum_required(VERSION 3.10)
project(ToyMonteCarlo CXX)

# 1. C++ Standard & Build Type
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimised build by default: validation suites run ~10^5-10^6 tracks
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ==============================================================================
# 2. ROOT Integration
# ==============================================================================
# Requires environment variables set by thisroot.sh.
# MathCore provides Landau density, CDF, and quantile functions (Ionization.cpp).
find_package(ROOT REQUIRED COMPONENTS Core Hist RIO Tree MathCore)
message(STATUS "Found ROOT Version: ${ROOT_VERSION}")

# ==============================================================================
# 3. Core Simulation Library (ToySimCore)
# ==============================================================================
# Compiled once from src/*.cpp and shared by the main demo and all validation apps.
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

add_library(ToySimCore STATIC ${SOURCES})
target_include_directories(ToySimCore PUBLIC include ${ROOT_INCLUDE_DIRS})
target_link_libraries(ToySimCore PUBLIC ${ROOT_LIBRARIES})

# Physics data directory (accessible from any working directory, like G4LEDATA):
#   - NIST XCOM photon attenuation tables (xcom_Z*.txt)
#   - ENDF/B-VIII.0 atomic relaxation data (atom-ZZZ_Sy_000.endf: MF=28)
#   - ENDF/B-VIII.0 photo-atomic data (photoat-ZZZ_Sy_000.endf: MF=23 subshells,
#     coherent/incoherent cross sections, and MF=27 F(x), S(x), f'(E), f''(E))
target_compile_definitions(ToySimCore PUBLIC TOYSIM_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# ==============================================================================
# 4. Main Simulation Application
# ==============================================================================
# Runs a single verbose event with full geometry, physics list, and ROOT output.
add_executable(SimulationApp app/main.cpp)
target_link_libraries(SimulationApp PRIVATE ToySimCore)

# ==============================================================================
# 5. Physics Validation Executables
# ==============================================================================

# 5a. Muon energy loss and range in lead vs PDG tables
add_executable(ValidationApp app/validateMuonLead.cpp)
target_link_libraries(ValidationApp PRIVATE ToySimCore)
target_compile_definitions(ValidationApp PRIVATE PDG_TABLE_PATH="${CMAKE_SOURCE_DIR}/data/muE_lead_Pb.txt")

# 5b. Electromagnetic physics & shower validation vs NIST XCOM, ESTAR, and PDG formulas
add_executable(ValidationEMApp app/validateEM.cpp)
target_link_libraries(ValidationEMApp PRIVATE ToySimCore)

# 5c. Particle decays: Michel spectrum, branching ratios, decay-time law
#     (in flight and at rest), kinematics, and neutrino energy balance
add_executable(ValidationDecayApp app/validateDecay.cpp)
target_link_libraries(ValidationDecayApp PRIVATE ToySimCore)

# 5d. Delta-ray production: differential cross sections vs numerical integration,
#     sampled spectra, kinematics, and cut-independence of mean dE/dx and range
add_executable(ValidationDeltaApp app/validateDelta.cpp)
target_link_libraries(ValidationDeltaApp PRIVATE ToySimCore)

# 5e. Atomic relaxation in lead: ENDF/B-VIII.0 binding energies and
#     K-line energies vs reference tables, subshell selection vs XCOM edge jumps,
#     cascade sampling, event-by-event energy conservation, and K X-ray escape
add_executable(ValidationRelaxApp app/validateRelaxation.cpp)
target_link_libraries(ValidationRelaxApp PRIVATE ToySimCore)

# 5f. Coherent (Rayleigh) & bound incoherent (Compton) scattering:
#     F(x)/S(x) limits and anomalous differential integrals vs ENDF MF=23 MT=502/504,
#     cross sections vs NIST XCOM, angular and energy-fraction sampling (quantile
#     chi2 and forward suppression), energy conservation, and thin-slab transmission
add_executable(ValidationScatteringApp app/validateRayleighScattering.cpp)
target_link_libraries(ValidationScatteringApp PRIVATE ToySimCore)

# 5g. Stopped negative muon capture & bound decay-in-orbit:
#     Suzuki et al. capture rates & Huff factors, Pb/C effective lifetimes &
#     capture fractions, mu+ non-capture, prompt K_alpha/L_alpha X-ray yields,
#     H 2-body kinematics, Pb neutron multiplicity, and energy/momentum conservation
add_executable(ValidationMuonCaptureApp app/validateMuonCapture.cpp)
target_link_libraries(ValidationMuonCaptureApp PRIVATE ToySimCore)

# 5h. [Stage 9] Landau-Pomeranchuk-Migdal (LPM) & dielectric suppression:
#     Migdal phi(s)/G(s)/xi(s) identities, SLAC E-146 25 GeV bremsstrahlung,
#     PeV pair suppression, quantile chi2 sampling, and 100 TeV shower elongation
add_executable(ValidationLPMApp app/validateLPM.cpp)
target_link_libraries(ValidationLPMApp PRIVATE ToySimCore)

# ==============================================================================
# 6. CTest Registration
# ==============================================================================
enable_testing()
add_test(NAME Stage7b_Relaxation COMMAND ValidationRelaxApp quick)
add_test(NAME Stage7c_Scattering COMMAND ValidationScatteringApp quick)
add_test(NAME Stage8_MuonCapture COMMAND ValidationMuonCaptureApp quick)
add_test(NAME Stage9_LPM COMMAND ValidationLPMApp quick)