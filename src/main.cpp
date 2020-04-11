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
#include "Potential.h"

#include "Test.h"

typedef int mode_t;

int main() {

	//Parameters testParameters = Parameters();
	//InitialConditions testInitialConditions = InitialConditions(&testParameters);
	//Potential testPotential = Potential(Vec3D(0, 0, 0));
	//double boxSize = 0.002;
	//std::vector<Star*> stars = testInitialConditions.initDiskStars(0, Vec3D(0.1, boxSize, 0), Vec3D(0.1+boxSize, 0, 0), 0.1, &testPotential);
	//InOut::write(stars, "testInitialConditions.dat");

	Test::initialConditionsSurfaceDensityBulge();

	return 0;

	Parameters parameters = Parameters();

	Database db = Database();
	db.open();
	db.setup();

	int selection;
	int simulationID = -1;
	std::cout << "[1] Load Simulation\n[2] New Simulation\n[3] Exit" << std::endl;
	std::cin >> selection;
	std::cin.clear();
	if (selection == 1) {
		std::vector<SimulationData> simulations = db.selectSimulations();
		std::cout << "Available Simulations:" << std::endl;
		for (SimulationData sim : simulations) {
			std::cout << sim.print();
		}
		std::cout << "Input ID to select simulation" << std::endl;
		std::cin >> selection;
		std::cin.clear();
		std::cout << "Running analysis on selected simulation" << std::endl;
		Analysis analysis = Analysis(parameters);
		int analysisID = db.insertAnalysis(selection, analysis);
		std::vector<int> timeSteps = db.selectTimesteps();
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
				db.insertAnalysisdtVelocity(analysisID,timeStep, analysis.average(velocities));
			}
		}
		std::cout << "Energy analysis done" << std::endl;
	}
	else if(selection == 2){
		Parameters* parameters = new Parameters();
		simulationID = db.insert(parameters);
		Simulation simulation = Simulation(simulationID);
		std::cout << "New simulation created: " << simulationID << std::endl;
		std::cout << "Starting simulation" << std::endl;
		simulation.run();
	}
	else if (selection == 3) {
		return 0;
	}

	

	//Vec3D test = Vec3D(2, 3, 4);
	//Vec3D normalVector = Vec3D(0, 0, 1);
	//Vec3D xy = Vec3D::projection(test,normalVector);

	//std::cout << xy.print() << std::endl;
	return 0;

	//Simulation
	//int nStars = 100;
	//double boxLength = 1; //[pc]
	//double dt = 0.01;
	//int nTimesteps = 100;
	//bool energyAnalysis = true;

	Analysis analysis = Analysis(parameters);


	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;

	
	std::cin.get();
	return 0;
}