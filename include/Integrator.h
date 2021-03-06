/**
 * Collection of Integrators. Calling one of the member functions will update coordinates and velocities of the stars by one timestep based on the given acceleration.
 *
 * @author Alarich Herzner
 * @version 0.10 15.05.2020
*/

#pragma once

#include <vector>

#include "Star.h"
#include "Node.h"
#include "MWPotential.h"

class Integrator
{
public:
    /** @brief Timestepsize used for integration. Dimension is day since the expected timeframe of interest is in range 0-10 years */
	double dt = 1; // [day]
    double dt2 = 0.5; //dt/2
    double dts; //[s]
    double dts2; //s/2
    double positionDetla; // dt*kmInpc*secInDay
    double positionDetla2; // 0.5*positionDetla
    /**@brief amount of seconds in one day. Needed because velocity is given in km/s and dt in days*/
    static double secInDay;
    /**@brief amount of pc in one km*/
    //static double kmInpc;

	Integrator(double dt = 0.01);

    /**
     @brief Implementation of the (in-)famous Euler algorithm (first-order method!). Performs one timestep per call. Updates the stars coordinates and velocity. 
     @param stars Vector of star pointers. All elements are updated
     @param dt Optional parameter timestepsize (only use this if timestepsize varies during simulation). If not passed, the timestepsize given on construction of the class instance is used.
     */
	void euler(std::vector<Star*> stars, double dt= 0);
    //https://math.stackexchange.com/questions/2023819/using-the-runge-kuttas-method-to-solve-a-2nd-derivative-question

    /**
    @brief Implementetion of "the" Runge Kutta algorithm.
    */
    void RK4(std::vector<Star*> stars, Node* root, MWPotential* potential);
    void RK4(std::vector<Star*> stars, MWPotential* potential);

    void Leapfrog(std::vector<Star*> stars, Node* root, MWPotential* potential);
    void Leapfrog(std::vector<Star*> stars, MWPotential* potential);
private:

};

