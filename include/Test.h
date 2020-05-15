#pragma once

#include <chrono>
#include <math.h>
#include <filesystem>

#include "Potential.h"
#include "InOut.h"
#include "Parameters.h"
#include "InitialConditions.h"
#include "Analysis.h"
#include "ProgressBar.h"
#include "Plot.h"

class Test : Parameters {
private: 
	static void pythonScript(std::string fileName);
	static std::string absolutePath;

public:

	static std::string path;

	static void sampleFieldStarPositions(int nStars = 5000);
	static void sampleFieldStarPositionsOutput(std::string path = "", int nStars = 5000);

	static void massDistribution(double z = 1000, double dx = 100);
	static void massDistributionDiskOutput(std::string path = "", double z = 1000, double dx = 100); // [z] = [pc]
	static void massDistributionBulgeOutput(std::string path = "", double z = 1000, double dx = 100); // [z] = [pc]

	static void potentialCircularVelocity();
	static void potentialCircularVelocityOutput(std::string path = "");

	static void velocityBulge();

	void initialConditionsMassSalpeterOutput(int nStars = 10000);
	static void initialConditionsMassBulgeOutput(double totalMass);
	static void potentialSurfaceDensityBulge();
	static void potentialSurfaceDensityDisk();
	static void initialConditionsSampleDisk(); //wip

	static void massDistributionTimer();

	static void velocityBulgeR();
	static void initialConditionsSampleBulgeVelocity();

	static void escapeVelocity();

	void initialConditionsInitFieldStars();

	static void velocityDispersionBulge();
};