/**
 * Barens Hut Simulation
 *
 * Stars are initialized (mass, possition- and velocity distribution). Octree is initialized based on the stars.
 * Stars accelerations are updated via the octree. Velocity and possition of the stars are updated via an integrator.
 * Parts of the code can be run in parallel via OpenMP
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#include <iostream>
#include <random>
#include <vector>
#include <chrono> //for timer
#include <direct.h> //_mkdir

#include "Star.h"
#include "Node.h"
#include "InOut.h"
#include "Integrator.h"
#include "InitialConditions.h"
#include "Analysis.h"
#include "Database.h"
#include "Simulation.h"
#include "MWPotential.h"
#include "WangPotential.h"

#include "Test.h"
#include "Plot.h"

typedef int mode_t;
//namespace plt = matplotlibcpp;

//global parameters
bool debug = false;

int main() {

	//InitialConditions testInitialConditions = InitialConditions();
	//MWPotential testPotential = MWPotential(Vec3D(0, 0, 0));
	//double boxSize = 0.002;
	//std::vector<Star*> stars = testInitialConditions.initDiskStars(0, Vec3D(0.1, boxSize, 0), Vec3D(0.1+boxSize, 0, 0), 0.1, &testPotential);
	//InOut::write(stars, "testInitialConditions.dat");

	//Test::potentialSurfaceDensityBulge();
	//Test::potentialSurfaceDensityDisk();

	//Test::potentialCircularVelocityOutput();
	//Test::velocityBulge();
	//Test::potentialCircularVelocityOutput();
	//Test::initialConditionsSampleBulgeVelocity();
	//Test::escapeVelocity();
	//Test::initialConditionsInitFieldStars();
	//MWPotential::generateVelocityDistributionBulgeLookupTable(25000);
	//std::vector<std::vector<double>> test = InOut::readDoubleMatrix("velocityDistributionBulgeTable.dat");

	Database db = Database();
	db.open();
	db.setup();

	while (true) {
		int selection;
		int simulationID = -1;
		std::cout << "[1] Load Simulation\n[2] New Simulation\n[3] Run Tests\n[4] Exit" << std::endl;
		std::cin >> selection;
		std::cin.clear();
		if (selection == 1) {
			if(!db.printSimulations())
				continue;
			std::cout << "Input ID to select simulation" << std::endl;
			int simulationID = 0;
			std::cin >> simulationID;
			std::cin.clear();
			db.selectSimulation(simulationID);
			Simulation simulation = Simulation(simulationID,&db);
			std::cout << "[1] Ouput\n[2] Analysis\n[3] Generate 2D" << std::endl;
			std::cin >> selection;
			std::cin.clear();
			if (selection == 1) {
				std::string directory = "Simulation" + std::to_string(simulation.getID());
				directory = InOut::makeDirectory(directory);
				std::cout << "Files will be written to: " << directory << std::endl;
				db.outputStars(simulation.getID(), directory + "/stars",false,true,true);
				std::cout << "done" << std::endl;
			}
			else if(selection==2) {
				std::cout << "Running analysis on selected simulation" << std::endl;
				Analysis analysis = Analysis();
				int analysisID = db.insertAnalysis(simulationID, analysis);
				/*std::vector<int> timeSteps = db.selectTimesteps();
				for (int timeStep : timeSteps) {
					std::vector<Star*> stars = db.selectStars(selection, timeStep);
					if (analysis.getbEnergy()) {
						db.insertAnalysisdtEnergy(analysisID, timeStep, analysis.kineticEnergy(stars), analysis.potentialEnergy(stars));
					}
					if (analysis.getbAverageVelocity()) {
						std::vector<Vec3D*> velocities = {};
						for (Star* star : stars) {
							velocities.push_back(&star->velocity);
						}
						double test = analysis.average(velocities);
						db.insertAnalysisdtVelocity(analysisID, timeStep, analysis.average(velocities));
					}
				}
				std::cout << "Energy analysis done" << std::endl;*/
				std::cout << "Running cluster analysis ..." << std::endl;
				analysis.cluster(db.selectPoints(simulationID,0,2));
				std::cout << "Cluster analysis done" << std::endl;
			}
			else {
				std::cout << "generating 2D positions and velocities ..." << std::endl;
				db.generate2D(simulation.getID());
				std::cout << "done\n" << std::endl;
			}
		}
		else if (selection == 2) {
			simulationID = db.insertSimulation();
			Simulation simulation = Simulation(simulationID, &db);
			std::cout << "New simulation created. ID = " << simulationID << std::endl;
			std::cout << "Starting simulation" << std::endl;
			simulation.run();
		}
		else if (selection == 3) {
			Test test = Test();
			Test::potentialCircularVelocity();
			test.transformation();
			Test::massDistribution(500,15000);
			Test::sampleFieldStarPositions(200);
			//Test::velocityBulgeR();
			Test::bulgeMass();
			test.velocityBulge();
			//Test::wangPositions();
			//Test::checkBrokenPowerLaw();
			std::cout << WangPotential::ANLM(1, 0, 0) << std::endl;
			std::cout << WangPotential::totalMass(-5e3, 5e3) << std::endl;
		}
		else if (selection == 4) {
			return 0;
		}
	}
	return 0;
}