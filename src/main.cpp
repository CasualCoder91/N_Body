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

#include "Star.h"
#include "Node.h"
#include "InOut.h"
#include "Integrator.h"
#include "InitialConditions.h"
#include "Analysis.h"
#include "Database.h"
#include "Simulation.h"

void runSimulation(Simulation& simulation);

int main() {

	Database db = Database();
	db.open("Data/test.db");
	db.setup();

	int selection;
	int simulationID = -1;
	std::cout << "[1] Load Simulation (coming soon)\n[2] New Simulation" << std::endl;
	std::cin >> selection;
	std::cin.clear();
	if (selection == 1) {
		std::vector<Simulation> simulations = db.selectSimulations();
		std::cout << "Available Simulations:" << std::endl;
		for (Simulation sim : simulations) {
			std::cout << sim.print();
		}
		std::cout << "Input ID to select simulation" << std::endl;
		std::cin >> selection;
		std::cin.clear();
		std::cout << "Running analysis on selected simulation" << std::endl;
		Analysis analysis = Analysis();
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
		Parameters parameters = Parameters();
		simulationID = db.insert(parameters);
		Simulation simulation = Simulation(simulationID);
		std::cout << "New simulation created: " << simulationID << std::endl;
		std::cout << "Starting simulation" << std::endl;
		runSimulation(simulation);
	}

	

	//Vec3D test = Vec3D(2, 3, 4);
	//Vec3D normalVector = Vec3D(0, 0, 1);
	//Vec3D xy = Vec3D::projection(test,normalVector);

	//std::cout << xy.print() << std::endl;
	return 0;

	//Simulation
	//int n_Stars = 100;
	//double boxLength = 1; //[pc]
	//double dt = 0.01;
	//int nTimesteps = 100;
	//bool energyAnalysis = true;

	Analysis analysis = Analysis();


	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;

	
	std::cin.get();
	return 0;
}

void runSimulation(Simulation& simulation){
	Database database = Database();
	int nextStarIndex = database.selectLastID("star")+1;
	//Init stars
	InitialConditions initialConditions = InitialConditions();
	std::vector<Star*> stars = initialConditions.initStars(nextStarIndex,simulation.getN_Stars());
	double totalMass = initialConditions.initialMass(stars);
	initialConditions.plummerSphere(stars, 1, totalMass);
	database.insertStars(simulation.getID(), stars, 0);

	//Integrate
	Integrator rk4 = Integrator(simulation.getdt());
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	for (int i = 1; i < simulation.getNTimesteps(); i++) {

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::findCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr);
		for (Star* star : stars) {
			root.insert(star);
		}

		root.calculateMassDistribution();

		rk4.euler(stars, &root);

		if (i % 10 == 0) {
			//InOut::writeWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::writeAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			database.timestep(i, stars);
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(stars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;

}