#pragma once

#include <string>
#include <vector>
#include "Configuration.h"

namespace Constants {

	static const std::string filePath = "./simulation.cfg";
	static Configuration config = Configuration(filePath);

	static const double softening = config.GetDouble(std::string("softening"));
	static const double precission = config.GetDouble(std::string("precission"));
	static const std::string title;
	static const double boxLength = config.GetDouble(std::string("boxLength"));
	static const double dt = config.GetDouble(std::string("dt")); //[day]
	static const int nTimesteps;
	static const int outputTimestep;
	static const int simulationID;

	//cluster mass
	static const double minMass;
	static const double maxMass;
	static const double alpha;
	static const std::vector<double> massLimits; //broken powerlaw
	static const std::vector<double> exponents; //broken powerlaw

	//View cone
	static const double angleOfView; //rad
	static const double dx; //pc
	static const double distance; //pc
	static const Vec3D viewPoint;
	static const Vec3D focus;

	//General Parameters, todo: put into separate class
	/** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
	static const double G;// = 4.483e-3;
	static const int nStars;
	static const Vec3D clusterLocation;
}