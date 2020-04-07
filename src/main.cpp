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

typedef int mode_t;

int main() {

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

//void runSimulation(Simulation& simulation,Parameters parameters){
//	Database database = Database();
//	int nextStarIndex = database.selectLastID("star")+1;
//	//Init stars
//	InitialConditions initialConditions = InitialConditions(parameters);
//	std::vector<Star*> stars = initialConditions.initStars(nextStarIndex,simulation.getNStars());
//	double totalMass = initialConditions.initialMass(stars);
//	initialConditions.plummerSphere(stars, 1, totalMass);
//	database.insertStars(simulation.getID(), stars, 0);
//
//	//Integrate
//	Integrator rk4 = Integrator(simulation.getdt());
//	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
//	for (int i = 1; i < simulation.getNTimesteps(); i++) {
//
//		Vec3D tlf = Vec3D(), brb = Vec3D();
//		Node::findCorners(tlf, brb, stars);
//		Node root = Node(tlf, brb, nullptr, parameters);
//		for (Star* star : stars) {
//			root.insert(star);
//		}
//
//		root.calculateMassDistribution();
//
//		rk4.euler(stars, &root);
//
//		if (i % 10 == 0) {
//			//InOut::writeWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
//			//InOut::writeAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
//			database.timestep(i, stars);
//		}
//	}
//	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
//	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
//	//InOut::write(stars,"stars.dat");
//	//InOut::write(&root);
//	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
//	std::cout << "done" << std::endl;
//
//}

void testfrequencyDistribution(){
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<Vec3D> Output;

	double z = 1; //kpc

	for (double x = -30; x < 30; x+=0.5) {
		std::cout << "x=" << x << std::endl;
		for (double y = -30; y < 30; y+=0.5) {
			double starMass = potential.frequencyDistribution(Vec3D(x, y, z));
			Output.push_back(Vec3D(x, y, starMass));
		}
	}

	InOut::write(Output, "frequencyDistribution_z" + std::to_string(z) + ".dat");

}

void samplePotential() {
	Potential potential = Potential(Vec3D(0, 0, 0));

	std::vector<Vec3D> positionsDisk;
	std::vector<Vec3D> positionsBulge;
	for (int i = 0; i < 10000; i++) {
		positionsDisk.push_back(potential.sampleDisk(-40, 40, -40, 40, -20, 20));
		positionsBulge.push_back(potential.sampleBuldge(-4, 4, -4, 4, -4, 4));
	}
	InOut::write(positionsDisk, "potentialTestDiskSample.dat");
	InOut::write(positionsBulge, "potentialTestBulgeSample.dat");
}

void potentialVelocityOutput() {
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<double> positions;
	std::vector<double> velocities;
	std::vector<double> disk;
	std::vector<double> blackHole;
	std::vector<double> buldge;
	std::vector<double> halo;

	for (double i = 0.1; i < 30; i += 0.1) {
		positions.push_back(i);
		Vec3D position = Vec3D(i, 0, 0);
		velocities.push_back(potential.circularVelocity(&position));
		disk.push_back(potential.circularVelocityDisk(&position));
		blackHole.push_back(potential.circularVelocityBlackHole(&position));
		buldge.push_back(potential.circularVelocityBulge(&position));
		halo.push_back(potential.circularVelocityHalo(&position));
	}

	InOut::write(positions, velocities, "potentialTest.dat");
	InOut::write(positions, disk, "potentialTestDisk.dat");
	InOut::write(positions, blackHole, "potentialTestBlackHole.dat");
	InOut::write(positions, buldge, "potentialTestBuldge.dat");
	InOut::write(positions, halo, "potentialTestHalo.dat");
}