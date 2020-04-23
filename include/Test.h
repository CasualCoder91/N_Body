#pragma once

#include <chrono>

#include "Potential.h"
#include "InOut.h"
#include "Parameters.h"
#include "InitialConditions.h"

class Test {
public:
	static void samplePotentialOutput(int nStars = 5000);
	void potentialCircularVelocityOutput();
	void testfrequencyDistribution();
	void initialConditionsMassSalpeterOutput(int nStars = 10000);
	static void initialConditionsMassBulgeOutput(double totalMass);
	static void potentialSurfaceDensityBulge();
	static void potentialSurfaceDensityDisk();
	static void initialConditionsSampleDisk(); //wip
	//Potential
	static void massDistributionDiskOutput(double z = 1); // [z] = [kpc]
	static void massDistributionBulgeOutput(double z = 1); // [z] = [kpc]


	//Timer
	static void massDistributionTimer();


};