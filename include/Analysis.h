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

class Analysis
{
private:
    bool bEnergy;
    bool bAverageVelocity;
    bool bAverage2DVelocity;
    double G;
public:

    Analysis(Parameters parameters);

    std::vector<double> totE;
    std::vector<double> potE;
    std::vector<double> kinE;
    std::vector<double> time;

    /**
     @brief The total potential energy of the given stars is calculated and stored into potE
     @param stars Vector of star pointers of which the potential energy is calculated.
     @return The potential energy of the given stars
     */
	double potentialEnergy(std::vector<Star*>& stars);
    /**
     @brief The total kinetic energy of the given stars is calculated and stored into kinE
     @param stars Vector of star pointers of which the kinetic energy is calculated.
     @return The kinetic energy of the given stars
     */
	double kineticEnergy(std::vector<Star*>& stars);
    /**
     @static
     @brief The scaling of the running time is tracked given a specific integrator depending on the number of stars and output into Output/NlogNtest.dat.
     @param maxNStars the maximum amount of stars to calculate the running time for.
     @param nTimesteps the amount of timesteps used for each number of stars.
     @param integrator the integrator used for integration over time.
     @param parameters Parameters read from cfg, passed as argument to avoid unnecessary reads. 
     */
	void static scaling(int maxNStars, int nTimesteps, Integrator& integrator, Parameters* parameters);

    bool getbEnergy();

    bool getbAverageVelocity();

    bool getbAverage2DVelocity();
    /**
     @static
     @brief calculates the average vector of the given \p vectors
     */
    static double average(std::vector<Vec3D*>& vectors);
    /**
     @static
     @brief calculates the dispersion (= variability/scatter/spread) of the given \p vectors
     */
    static double dispersion(std::vector<Vec3D*>& vectors);
    /** @brief saves the calculated energy values for each timestep to .dat Files (TotalEnergy.dat, KinetikEnergy.dat,PotentialEnergy.dat)*/
    void write();
};

