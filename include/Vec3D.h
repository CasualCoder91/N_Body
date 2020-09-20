/**
 * Implementation of a 3D Vector.
 * 
 * Position, velocity and acceleration of stars are stored as Vec3D.
 * All neccesary operations are supported via methods/operator overloading.
 *
 * @attention All member variables are public.
 * @author Alarich Herzner
 * @version 0.9 02.03.2020
*/

#pragma once
#include <cstdlib> //abs()
#include <string>
#include <random>

#include "Matrix.h"

class Vec3D{
public:
    /** @brief x value of the vector */
    double x;
    /** @brief y value of the vector */
    double y;
    /** @brief z value of the vector */
    double z;

    Vec3D();
    Vec3D(double x, double y, double z);
    ~Vec3D();
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
    static Vec3D randomAngles(double length);
    /**
     @brief Calculates lenght (euclidean norm).
     @return lenght of the vector.
     */
    double length();

    double phi();

    double theta();

    Vec3D cartesianToCylindrical();

    Vec3D cartesianToCylindricalV(Vec3D cartesianPos);

    Vec3D cartesianToSphericalV(Vec3D cartesianPos);

    /**
     @brief Calculates normalized vector with direction of caller.
     @return Vector with lenght=1 and direction of caller.
     */
    Vec3D normalize();
    /**
     @static
     @brief Calculates the cross product of two vectors.
     @param a, b Vectors of which the cross product is calculated.
     @return cross product.
     */
    static Vec3D crossProduct(const Vec3D* a, const Vec3D* b);
    /**
     @brief creates string containing coordinates of the vector. Used for output.
     @return coordinates of vector as sting in format: x,y,z.
     */
    std::string print();
    /**
     @static
     @brief Calculates distance (euclidean norm) between given vectors.
     @param a, b Vectors of which the cross product is calculated.
     @return distance between given vectors.
     */
    static double distance(const Vec3D* a, const Vec3D* b);

    static double distance2(const Vec3D* a, const Vec3D* b);
    /**
    @brief Sets all member variables to 0.
    */
    void reset();
    /**
     @brief Operator overloading for division of a vector by a double value.
     */
    friend Vec3D operator / (const Vec3D& lhs, double& rhs);
    /**
    @brief Operator overloading for multiplication of a vector by a double value.
    @see Vec3D operator * (double& lhs, const Vec3D& rhs)
    */
    friend Vec3D operator * (const Vec3D& lhs, const double& rhs);
    /**
    @brief Operator overloading for multiplication of a double value by a vector.
    @see Vec3D operator * (const Vec3D& lhs, double& rhs)
    */
    friend Vec3D operator * (const double& lhs, const Vec3D& rhs);
    /**
    @brief Operator overloading for multiplication of a vector with a vector.
    */
    double operator * (const Vec3D& rhs);
    /**
    @brief Operator overloading for addition of two vectors.
    @see friend Vec3D operator + (Vec3D lhs, Vec3D const& rhs);
    */
    Vec3D& operator += (const Vec3D& rhs);
    /**
    @brief Operator overloading for addition of two vectors.
    @see Vec3D& operator += (const Vec3D& rhs);
    */
    friend Vec3D operator + (Vec3D lhs, Vec3D const& rhs);
    friend Vec3D operator - (Vec3D lhs, Vec3D const& rhs);
    static Vec3D projectPosition(const Vec3D& target, const Vec3D& lookAt, const Vec3D& origin, const double fovAngle);
    static Vec3D project(Vec3D& position, Vec3D& velocity,  Vec3D& lookAt, Vec3D& origin);
    static Vec3D projectVelocity(const Vec3D& target, const Vec3D& lookAt, const double fovAngle);
};

