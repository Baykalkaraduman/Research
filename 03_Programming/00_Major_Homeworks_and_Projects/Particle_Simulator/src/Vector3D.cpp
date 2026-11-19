#include "Vector3D.h"

#include <cmath>

// ==============================================================================
// rotateUz
// ==============================================================================
// In the local frame the old direction is (0, 0, 1). A vector given there, e.g.
// (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta)), is turned so that the
// local z axis lands on newUz. The rotation axis is perpendicular to both z and newUz.
Vector3D& Vector3D::rotateUz(const Vector3D& newUz) {
    const double u1 = newUz.x;
    const double u2 = newUz.y;
    const double u3 = newUz.z;
    double up = u1 * u1 + u2 * u2;

    if (up > 0.0) {
        up = std::sqrt(up);
        const double px = x, py = y, pz = z;
        x = (u1 * u3 * px - u2 * py) / up + u1 * pz;
        y = (u2 * u3 * px + u1 * py) / up + u2 * pz;
        z = -up * px + u3 * pz;
    } else if (u3 < 0.0) {
        // newUz = (0, 0, -1): a rotation by pi around the y axis
        x = -x;
        z = -z;
    }
    // newUz = (0, 0, +1): local and global frames are the same, nothing to do
    return *this;
}
