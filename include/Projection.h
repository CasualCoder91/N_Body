#pragma once

#include "Vec3D.h"
#include "Matrix.h"

class Projection {
    //static double kmInpc;
public:
    static Vec3D projectPosition(const Vec3D& target, const Vec3D& lookAt, const Vec3D& origin, const double fovAngle);
    static void project(Vec3D& position, Vec3D& velocity, Vec3D& lookAt, Vec3D& origin);
    static Vec3D projectVelocity(const Vec3D& target, const Vec3D& lookAt, const double fovAngle);
};