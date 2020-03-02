/**
 * Implementation of a 3D Vector.
 * 
 * Position, velocity and acceleration of stars are stored as Vec3D.
 * All neccesary operations are supported via methods/operator overloading.
 *
 * @author Alarich Herzner
 * @version 0.9 02.03.2020
*/

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
    /** 
    @static
    @brief Calculates vector halve way between input vectors @p a and @p b. 
    Used to get the center point of octree nodes: #Node::center.
    @param a, b the vectors of which the middle is evaluated
    */
    static Vec3D center(Vec3D a, Vec3D b);
    /**
    @static
    @brief Calculates 3D vector pointing in random direction with given lenght (euclidean norm).
    @param length Euclidean norm of the returned vector.
    @return the calculated 3D vector.
    @see randomAngles()
    */
    static Vec3D randomVector(double length);
    /**
     @static
     @brief Calculates 3D vector pointing in random direction with given lenght (euclidean norm).
     @param length Euclidean norm of the returned vector.
     @return the calculated 3D vector.
     @see randomAngles()
     */
    static Vec3D randomAngles(double r);
    double length();
    Vec3D normalize();
    static Vec3D crossProduct(Vec3D v1, Vec3D v2);
    std::string print();
    static double distance(Vec3D* a, Vec3D* b);
    void reset();
    Vec3D operator / (double& rhs);
    Vec3D operator * (double& rhs);
    friend Vec3D operator*(float f, const Vec3D& v);
    double operator * (Vec3D& rhs);
    Vec3D operator + (Vec3D const& rhs);
    Vec3D& operator += (const Vec3D &rhs);
};

