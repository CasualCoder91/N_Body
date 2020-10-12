#pragma once

#include <string>
#include <vector>
#include "Configuration.h"

namespace Constants {

	extern std::string filePath;
	extern Configuration config;

	extern double softening;
	extern double precission;
	extern std::string title;
	extern double boxLength;
	extern double dt; //[day]
	extern int nTimesteps;
	extern int outputTimestep;
	extern int simulationID;

	//cluster mass
	extern double minMass;
	extern double maxMass;
	extern double alpha;
	extern std::vector<double> massLimits; //broken powerlaw
	extern std::vector<double> exponents; //broken powerlaw

	//View cone
	extern double angleOfView; //rad
	extern double dx; //pc
	extern double distance; //pc
	extern Vec3D viewPoint;
	extern Vec3D focus;

	//General Parameters, todo: put into separate class
	/** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
	extern double G;// = 4.483e-3;
	extern int nStars;
	extern Vec3D clusterLocation;

	//Analysis Parameters
	extern bool bEnergy;
	extern bool bAverageVelocity;
	extern bool bAverage2DVelocity;
}