#pragma once

#include <random>

class RandomEngine {
public:
    // std::mt19937_64::result_type is named as seedType for convenience
    using seedType = std::mt19937_64::result_type;

private:

// Defining the default seed
static constexpr seedType defaultSeed = 1722;

// Random number generation engine:
// Initialized with "defaultSeed"
static inline std::mt19937_64 engine{defaultSeed};

// Seed is saved as "theSeed" and it is initialized with the "defaultSeed"
static inline seedType theSeed = defaultSeed;
    
public:

    // ===================================================
    // GETTERS
    // ===================================================

    // Gets the seed
    static seedType getSeed() {return theSeed;} 

    // ===================================================
    // SETTERS
    // ===================================================

    // Defines a new seed and gives it to the engine
    static void setSeed(seedType newSeed) {theSeed = newSeed; engine.seed(newSeed);}

    // ===================================================
    // UTILITY
    // ===================================================
    // Uniform number generation
    static double uniNumGen() {return std::uniform_real_distribution<double>(0.0, 1.0)(engine);}

    // Gaussian random number
    static double randGauss(double mean, double sigma) {return std::normal_distribution<double>(mean, sigma)(engine);}

    // Gamma-distributed random number
    static double randGamma(double shape, double scale) {return std::gamma_distribution<double>(shape, scale)(engine);}

    // Direct engine access, for other std:: distributions (normal, exponential, ...)
    static std::mt19937_64& getEngine() {return engine;}

};