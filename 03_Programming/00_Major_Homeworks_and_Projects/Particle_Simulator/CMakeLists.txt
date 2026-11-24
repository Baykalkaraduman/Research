# 1. Define minimum CMake version and Project Name
cmake_minimum_required(VERSION 3.10)
project(ToyMonteCarlo CXX)

# 2. Modern C++ Standard (Required for smart pointers and ROOT 6+)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optimised build by default: the validation runs ~10^5 tracks
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ==============================================================================
# 3. ROOT INTEGRATION
# ==============================================================================
# Find ROOT packages using the environment variables set by thisroot.sh
# + MathCore: the Landau density, CDF and quantile (Ionization.cpp)
find_package(ROOT REQUIRED COMPONENTS Core Hist RIO Tree MathCore)

# Print a message to the terminal to confirm ROOT was found
message(STATUS "Found ROOT Version: ${ROOT_VERSION}")

# 4. Gather Source Files (This automatically picks up AnalysisManager.cpp)
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

# ==============================================================================
# 5. The simulation code as a library, used by two programs
# ==============================================================================
# Before: one executable compiled from app/main.cpp + src/*.cpp. Now there are two
# programs (the demo and the validation); a static library compiles src/ once.
add_library(ToySimCore STATIC ${SOURCES})
target_include_directories(ToySimCore PUBLIC include ${ROOT_INCLUDE_DIRS})
target_link_libraries(ToySimCore PUBLIC ${ROOT_LIBRARIES})

# 6. The demo (one verbose event)
add_executable(SimulationApp app/main.cpp)
target_link_libraries(SimulationApp PRIVATE ToySimCore)

# 7. Validation against the PDG muon table for lead
add_executable(ValidationApp app/validateMuonLead.cpp)
target_link_libraries(ValidationApp PRIVATE ToySimCore)
# Path of the PDG table, so the program finds it from any working directory
target_compile_definitions(ValidationApp PRIVATE PDG_TABLE_PATH="${CMAKE_SOURCE_DIR}/data/muE_lead_Pb.txt")

# 8. Data tables used by the physics itself, found from any working directory
#    (like Geant4's G4LEDATA): XCOM photon cross sections and, since stage 7b,
#    the ENDF/B-VIII.0 files for lead (atom-082_Pb_000.endf: atomic relaxation;
#    photoat-082_Pb_000.endf: subshell photoelectric cross sections)
target_compile_definitions(ToySimCore PUBLIC TOYSIM_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# 9. Validation of the photon / electron physics and of the shower
#    against NIST XCOM, NIST ESTAR and the PDG shower formulas
add_executable(ValidationEMApp app/validateEM.cpp)
target_link_libraries(ValidationEMApp PRIVATE ToySimCore)

# 10. Validation of the decays: Michel spectrum, branching ratios,
#     kinematics, the decay-time law (in flight and at rest) and the energy
#     balance with neutrinos
add_executable(ValidationDecayApp app/validateDecay.cpp)
target_link_libraries(ValidationDecayApp PRIVATE ToySimCore)

# 11. Validation of the delta rays: cross sections vs numerical
#     integration, sampled spectra, collision kinematics, and that the mean loss
#     and the mean range do not depend on the production cut
add_executable(ValidationDeltaApp app/validateDelta.cpp)
target_link_libraries(ValidationDeltaApp PRIVATE ToySimCore)

# 12. [Stage 7b] Validation of the atomic relaxation: shell data vs reference
#     values, subshell selection vs the XCOM K-edge jump, the cascade sampling,
#     energy conservation and the K X-ray escape from a lead block vs the
#     analytic formula
add_executable(ValidationRelaxApp app/validateRelaxation.cpp)
target_link_libraries(ValidationRelaxApp PRIVATE ToySimCore)
