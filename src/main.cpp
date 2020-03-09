/**
 * Barnes–Hut simulation.
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

int main() {

	Database db = Database();
	db.open("Data/test.db");
	db.setup();

	int selection;
	int simulationID = -1;
	Parameters parameters = Parameters();
	std::cout << "[1] Load Simulation (coming soon)\n[2] New Simulation" << std::endl;
	std::cin >> selection;
	if (selection == 2) {
		simulationID = db.insert(parameters);
		std::cout << "New Simulation created: " << simulationID << std::endl;
	}

	//Init stars
	std::vector<Star*> stars = {};
	double totalMass = InitialConditions::initialMass(stars, parameters.getN_Stars());
	InitialConditions::plummerSphere(stars, 1, totalMass);
	db.insert(simulationID, stars);
	std::cin.get();
	//Vec3D test = Vec3D(2, 3, 4);
	//Vec3D normalVector = Vec3D(0, 0, 1);
	//Vec3D xy = Vec3D::projection(test,normalVector);

	//std::cout << xy.print() << std::endl;
	return 0;

	//Parameters
	//int n_Stars = 100;
	//double boxLength = 1; //[pc]
	//double dt = 0.01;
	//int nTimesteps = 100;
	//bool energyAnalysis = true;

	Analysis analysis = Analysis(parameters.getEnergyAnalysis());


	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;

	//Integrate
	Integrator rk4 = Integrator(parameters.getdt());

	//Analysis::scaling(5, 5, euler);
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	for (int i = 0; i < parameters.getNTimesteps(); i++) {

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::findCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr);
		for (Star* star : stars) {
			root.insert(star);
		}
		root.calculateMassDistribution();
		//#pragma omp parallel for //1:10
		//for (int i = 0; i < stars.size();++i){
		//	stars.at(i)->acceleration.reset();
		//	root.applyForce(stars.at(i));
		//}
		rk4.euler(stars,&root,parameters.getdt());

		if (i % 10 == 0) {
			//InOut::writeWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::writeAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			analysis.time.push_back(i);
			if (analysis.getbEnergy()) {
				double potentialEnergy = analysis.potentialEnergy(stars);
				double kineticEnergy = analysis.kineticEnergy(stars);
				std::cout << "stepnr: " << i << std::endl;
				std::cout << "Kinetic Energy: " + std::to_string(kineticEnergy) << std::endl;
				std::cout << "Potential Energy: " + std::to_string(potentialEnergy) << std::endl;
				std::cout << "Total Energy: " + std::to_string(kineticEnergy + potentialEnergy) << std::endl << std::endl;
			}
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(stars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
	std::cin.get();
	return 0;
}