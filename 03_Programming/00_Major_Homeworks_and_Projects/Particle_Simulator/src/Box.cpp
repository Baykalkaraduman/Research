#include "Box.h"
#include "GeoDefs.h"

#include <algorithm> // std::min, std::max, std::swap

// ==============================================================================
// Distance from an INSIDE point to the surface
// ==============================================================================
double Box::distanceToOut(const Vector3D& p, const Vector3D& v) const {
    const double pos[3]  = {p.x, p.y, p.z};
    const double dir[3]  = {v.x, v.y, v.z};
    const double half[3] = {hx, hy, hz};

    double distance = GeoDefs::kInfinity;
    for (int i = 0; i < 3; ++i) {
        // Moving towards +half or -half on this axis? (dir == 0: never leaves on this axis)
        if (dir[i] > 0.0) {
            distance = std::min(distance, (half[i] - pos[i]) / dir[i]);
        } else if (dir[i] < 0.0) {
            distance = std::min(distance, (-half[i] - pos[i]) / dir[i]);
        }
    }
    return std::max(distance, 0.0); // a point within tolerance outside gives a tiny negative value
}

// ==============================================================================
// Distance from an OUTSIDE point to the surface (kInfinity if it misses)
// ==============================================================================
double Box::distanceToIn(const Vector3D& p, const Vector3D& v) const {
    const double pos[3]  = {p.x, p.y, p.z};
    const double dir[3]  = {v.x, v.y, v.z};
    const double half[3] = {hx, hy, hz};

    double tNear = -GeoDefs::kInfinity; // the ray is inside ALL slabs between tNear ...
    double tFar  =  GeoDefs::kInfinity; // ... and tFar

    for (int i = 0; i < 3; ++i) {
        if (dir[i] == 0.0) {
            // Parallel to this slab: either always inside it, or never
            if (std::abs(pos[i]) > half[i]) {
                return GeoDefs::kInfinity;
            }
            continue;
        }
        double t1 = (-half[i] - pos[i]) / dir[i];
        double t2 = ( half[i] - pos[i]) / dir[i];
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        tNear = std::max(tNear, t1);
        tFar  = std::min(tFar, t2);
    }

    // No overlap of the three intervals, or the box is behind us
    if (tNear > tFar || tFar < 0.0) {
        return GeoDefs::kInfinity;
    }
    return std::max(tNear, 0.0);
}
