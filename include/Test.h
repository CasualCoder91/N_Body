#pragma once

#include <chrono>
#include <math.h>
#include <filesystem>

#include "MWPotential.h"
#include "InOut.h"
#include "InitialConditions.h"
#include "Analysis.h"
#include "ProgressBar.h"
#include "Plot.h"
#include "WangPotential.h"
#include "Matrix.h"
#include "Projection.h"

class Test {
private: 
	static void pythonScript(std::string fileName);
	static std::string absolutePath;
	MWPotential potential = MWPotential();
	InitialConditions initialConditions = InitialConditions(&potential);

public:
	Test();

	static std::string path;

	static void sampleFieldStarPositions(int nStars = 5000);
	static void sampleFieldStarPositionsOutput(std::string path = "", int nStars = 5000);

	static void massDistribution(double z = 1000, double dx = 100);
	static void massDistributionDiskOutput(std::string path = "", double z = 1000, double dx = 100); // [z] = [pc]
	static void massDistributionBulgeOutput(std::string path = "", double z = 1000, double dx = 100); // [z] = [pc]

	static void potentialCircularVelocity();
	static void potentialCircularVelocityOutput(std::string path = "");

	void initialConditionsMassSalpeterOutput(int nStars = 10000);
	static void initialConditionsMassBulgeOutput(double totalMass);
	static void potentialSurfaceDensityBulge();
	static void potentialSurfaceDensityDisk();
	static void initialConditionsSampleDisk(); //wip

	static void massDistributionTimer();

	static void velocityBulgeR();
	static void velocityDispersionBulgerGC();
	static void initialConditionsSampleBulgeVelocity();
	void velocityBulge();

	void velocityDisk();

	static void escapeVelocity();

	void initialConditionsInitFieldStars();
	std::vector<Star*> initBulgeStars(int& starID, Vec3D focus, Vec3D viewPoint, double distance, double angleOfView=0.005);

	static void bulgeMass();

	static void checkBrokenPowerLaw();

	static void transformation();

};