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
	extern double plummer_radius;
	extern double dt; //[day]
	extern int nTimesteps;
	extern int outputTimestep;

	//cluster
	extern std::vector<double> massLimits; //broken powerlaw
	extern std::vector<double> exponents; //broken powerlaw
	extern Vec3D clusterLocation;
	extern bool bMcLuster;
	extern std::string mcluster_filepath;

	//View cone
	extern double angleOfView; //rad
	extern double distance; //pc
	extern Vec3D viewPoint;
	extern Vec3D focus;

	//General Parameters
	/** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
	extern double G;// = 4.483e-3;
	extern double degInRad;
	extern double radInArcsec;
	extern double radmyrInArcsecyr;
	/**@brief 1km divided by 1pc*/
	extern double kmInpc;
	extern int nStars;

	//Analysis Parameters (deprecated)
	//extern bool bEnergy;
	//extern bool bAverageVelocity;
	//extern bool bAverage2DVelocity;

	//Transformation
	extern Vec3D positionSun; //in cylinder coordinates [pc,rad,pc]
	extern Vec3D velocitySun; //kms
	extern double circularVelocitySun; //kms

	//Math
	extern double pi;
	extern double pi2;
	extern double pi4;

	//Clustering
	extern double eps_magnitude; //[%] maximum change in magnitude to be considered the same star
}