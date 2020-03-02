#pragma once
#include <cstdlib> //abs()
#include <string>
#include <random>

class Vec3D{
public:
    double x;
    double y;
    double z;

    Vec3D();
    Vec3D(double x, double y, double z);
    static Vec3D Center(Vec3D a, Vec3D b);
    static Vec3D RandomVector(double length);
    static Vec3D RandomAngles(double r);
    double Length();
    Vec3D Normalize();
    static Vec3D Cross(Vec3D v1, Vec3D v2);
    std::string Print();
    static double Distance(Vec3D* a, Vec3D* b);
    void Reset();
    Vec3D operator / (double& rhs);
    Vec3D operator * (double& rhs);
    friend Vec3D operator*(float f, const Vec3D& v);
    double operator * (Vec3D& rhs);
    Vec3D operator + (Vec3D const& rhs);
    Vec3D& operator += (const Vec3D &rhs);
};

