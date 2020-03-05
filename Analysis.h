/**
 * All types of analysis are defined here.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once

#include <vector>
#include <chrono> //for timer

#include "Star.h"
#include "Parameters.h"
#include "Integrator.h"
#include "Node.h"
#include "InOut.h"
#include "InitialConditions.h"

class Analysis : Parameters
{
public:
    /**
     @static
     @brief The total potential energy of the given stars is calculated.
     @param stars Vector of star pointers of which the potential energy is calculated.
     @return The potential energy of the given stars
     */
	double static potentialEnergy(std::vector<Star*>& stars);
    /**
     @static
     @brief The total kinetic energy of the given stars is calculated.
     @param stars Vector of star pointers of which the kinetic energy is calculated.
     @return The kinetic energy of the given stars
     */
	double static kineticEnergy(std::vector<Star*>& stars);
    /**
     @static
     @brief The scaling of the running time is tracked given a specific integrator depending on the number of stars and output into Output/NlogNtest.dat.
     @param maxNStars the maximum amount of stars to calculate the running time for.
     @param nTimesteps the amount of timesteps used for each number of stars.
     @param integrator the integrator used for integration over time.
     */
	void static scaling(int maxNStars, int nTimesteps, Integrator& integrator);
};

