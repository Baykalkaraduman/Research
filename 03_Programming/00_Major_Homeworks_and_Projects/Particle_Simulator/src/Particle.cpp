#include "Particle.h"
#include <iomanip> // For std::setprecision

// ==============================================================================
// Constructor
// Using Member Initializer Lists for optimal performance and std::move for strings.
// ==============================================================================
Particle::Particle(int p_pdg_id, std::string p_name, double p_mass, double p_charge, double p_spin, double p_lifetime)
    : pdg_id(p_pdg_id), name(std::move(p_name)), mass(p_mass), charge(p_charge), spin(p_spin), lifetime(p_lifetime) {}

// ==============================================================================
// Utility Implementation
// ==============================================================================
bool Particle::isStable() const {
    // For convention, a negative lifetime implies stability.
    return lifetime < 0.0; 
}

// ==============================================================================
// Virtual Display Implementation
// ==============================================================================
void Particle::display() const {
    std::cout << "\n+--- Particle Definition: " << name << " -------------------+\n";
    std::cout << "| PDG ID       : " << pdg_id << "\n";
    std::cout << "| Mass         : " << std::fixed << std::setprecision(6) << mass << " GeV/c^2\n";
    std::cout << "| Charge       : " << std::defaultfloat << charge << " e\n";
    std::cout << "| Spin         : " << spin << "\n";
    
    if (isStable()) {
        std::cout << "| Lifetime     : Stable\n";
    } else {
        std::cout << "| Lifetime     : " << lifetime << " ns\n";
    }
    std::cout << "+---------------------------------------------------+\n";
}