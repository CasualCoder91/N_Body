/**
 * Collection of Integrators. Calling one of the member functions will update coordinates and velocities of the stars by one timestep based on the given acceleration.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once
#include "Star.h"
#include <vector>

class Integrator
{
public:
    /** @brief Timestepsize used for integration. Dimension is pc*s/km to minimize cutoff error caused by big differences in orders of magnetude */
	double dt = 0.1; //timestepsize

	Integrator(double dt = 0.01);

    /**
     @brief Implementation of the (in-)famous Euler algorithm (first-order method!). Performs one timestep per call. Updates the stars coordinates and velocity. 
     @param stars Vector of star pointers. All elements are updated
     @param dt Optional parameter timestepsize (only use this if timestepsize varies during simulation). If not passed, the timestepsize given on construction of the class instance is used.
     */
	void euler(std::vector<Star*> stars,double dt= 0); // one timestep with euler algorithm
};

