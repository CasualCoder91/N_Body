#pragma once

#include <cstdlib> //abs(),sqrt()

class Vec2D {
public:
    /** @brief x value of the vector */
    double x;
    /** @brief y value of the vector */
    double y;

    Vec2D();
    Vec2D(double x, double y);
    ~Vec2D();

    /**
     @brief Calculates lenght (euclidean norm).
     @return lenght of the vector.
     */
    double length() const;

};