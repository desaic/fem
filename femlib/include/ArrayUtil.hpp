#ifndef ARRAYUTIL_HPP
#define ARRAYUTIL_HPP
#include <vector>
#include "vecmath.h"

void BBox(const std::vector<Vector3f >& v,
    Vector3f & mn, Vector3f & mx);
#endif // ARRAYUTIL_HPP
