/**
 * All types of analysis are defined here.
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#pragma once

#include <vector>
#include <chrono> //for timer

#include <mlpack/core.hpp>
#include <mlpack/methods/dbscan/dbscan.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/methods/range_search/range_search.hpp>

#include "Star.h"
#include "Integrator.h"
#include "Node.h"
#include "InOut.h"
#include "InitialConditions.h"
#include "Vec2D.h"
#include <Point.h>
#include <Database.h>
#include <Vec3D.h>



class Analysis
{
private:

    static std::string path;
    int id;
    Database* database;

public: //variables

public: //methods

    Analysis(int id, Database* database);

    std::vector<double> totE;
    std::vector<double> potE;
    std::vector<double> kinE;
    std::vector<double> time;

    /**
     @brief The total potential energy of the given stars is calculated and stored into potE
     @param stars Vector of star pointers of which the potential energy is calculated.
     @return The potential energy of the given stars
     */
	double potentialEnergy(const std::vector<Star>& stars);
    /**
     @brief The total kinetic energy of the given stars is calculated and stored into kinE
     @param stars Vector of star pointers of which the kinetic energy is calculated.
     @return The kinetic energy of the given stars
     */
	double kineticEnergy(const std::vector<Star>& stars);
    /**
     @static
     @brief The scaling of the running time is tracked given a specific integrator depending on the number of stars and output into Output/NlogNtest.dat.
     @param maxNStars the maximum amount of stars to calculate the running time for.
     @param nTimesteps the amount of timesteps used for each number of stars.
     @param integrator the integrator used for integration over time.
     */
	void static scaling(int maxNStars, int nTimesteps, Integrator& integrator);

    void energy();

    /**
     @static
     @brief calculates the average vector of the given \p vectors
     */
    static double average(std::vector<Vec3D>& vectors);
    static double average(std::vector<Vec2D>& vectors);
    static double average(std::vector<double>& values);
    /**
     @static
     @brief calculates the average velocity of the given \p points
     */
    static double average(std::vector<Point>& points);

    /**
     @static
     @brief calculates the dispersion (= variability/scatter/spread) of the given \p vectors
     */
    static double dispersion(std::vector<Vec3D>& vectors, double average = -1);
    static double dispersion(std::vector<Vec2D>& vectors, double average = -1);
    static double dispersion(std::vector<double>& values, double average = -1);
    /**
     @static
     @brief calculates the velocity dispersion of the given \p points
     */
    static double dispersion(std::vector<Point>& points, double average = -1);

    /** @brief saves the calculated energy values for each timestep to .dat Files (TotalEnergy.dat, KinetikEnergy.dat,PotentialEnergy.dat)*/
    void write();


    /**
     @brief Generates HTP (2D) velocities for stars stored in the database. 
     @param observed set true to use stars observed with photutils, false to use simulated stars
     @param nTimesteps the amount of timesteps used for each number of stars.
     @param integrator the integrator used for integration over time.
     */
    void generateHTPVelocity(bool observed = false, bool force_correct_selection = false);

    void cluster(std::vector<Point>& points);

    /**
     @brief maps observed stars to simulated stars by setting star.fkStar of the observed to star.id to the simulated star.
     */
    void map_observed();
};

