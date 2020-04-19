/**
 * Collection of Integrators. Calling one of the member functions will update coordinates and velocities of the stars by one timestep based on the given acceleration.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once

#include <vector>

#include "Star.h"
#include "Node.h"

class Integrator
{
public:
    /** @brief Timestepsize used for integration. Dimension is pc*s/km to minimize cutoff error caused by big differences in orders of magnetude */
	double dt = 0.1; // [day]
    double dt2 = 0.05; //dt/2
    static double dayInSec;
    static double kmInpc;

    Vec3D* C[4]; //RK4
    Vec3D* K[4]; //RK4

	Integrator(double dt = 0.01);

    /**
     @brief Implementation of the (in-)famous Euler algorithm (first-order method!). Performs one timestep per call. Updates the stars coordinates and velocity. 
     @param stars Vector of star pointers. All elements are updated
     @param dt Optional parameter timestepsize (only use this if timestepsize varies during simulation). If not passed, the timestepsize given on construction of the class instance is used.
     */
	void euler(std::vector<Star*> stars, Node* root, double dt= 0);
    //https://math.stackexchange.com/questions/2023819/using-the-runge-kuttas-method-to-solve-a-2nd-derivative-question
    void RK4(std::vector<Star*> stars, Node* root, double dt = 0);
};

