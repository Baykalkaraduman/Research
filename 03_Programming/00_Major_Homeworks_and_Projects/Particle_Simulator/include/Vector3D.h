#pragma once

#include <iostream>
#include <cmath>

// ==============================================================================
// Struct: Vector3D
// Updated with basic kinematics math (magnitudes and spherical conversions).
// ==============================================================================
struct Vector3D {
    double x, y, z;

    Vector3D(double valX = 0.0, double valY = 0.0, double valZ = 0.0) 
        : x(valX), y(valY), z(valZ) {}

    // Calculate the magnitude (length) of the vector
    double mag() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // Generate a new vector from spherical coordinates
    // r: magnitude (momentum), theta: polar angle [0, pi], phi: azimuthal angle [0, 2pi]
    static Vector3D fromSpherical(double r, double theta, double phi) {
        return Vector3D(
            r * std::sin(theta) * std::cos(phi),
            r * std::sin(theta) * std::sin(phi),
            r * std::cos(theta)
        );
    }

    // ==============================================================================
    // Vector Arithmetics
    // ==============================================================================
    
    // Dot Product
    double dot(const Vector3D& v) const {return x*v.x + y*v.y + z*v.z;}

    // Cross Product
    Vector3D cross(const Vector3D& v) const {return Vector3D(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);}

    // Finding Unit Vector
    Vector3D unit() const {
        double m = mag(); // defining m so we don't need to calculate mag() twice
        if (m == 0) 
            return Vector3D(0.0, 0.0, 0.0);
        return *this/m;
    }

    // ------------------------------------------------------------------------------
    // Operator Overloading for calculations with vectors
    // ------------------------------------------------------------------------------

    // Multiplying Vectors with Scalars
    Vector3D operator*(double s) const {return Vector3D(x*s, y*s, z*s);}
    friend Vector3D operator*(double s, const Vector3D& v) {return v*s;} // now we can multiply in any order

    // Dividing Vectors with Scalars
    Vector3D operator/(double s) const {return Vector3D(x/s, y/s, z/s);}

    // Adding Vectors to each other
    Vector3D operator+(const Vector3D& v) const {return Vector3D(x+v.x, y+v.y, z+v.z);}

    // Subtracting Vectors from each other
    Vector3D operator-(const Vector3D& v) const {return Vector3D(x-v.x, y-v.y, z-v.z);}

    void display() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")";
    }
};