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

#include  "Star.h"
#include "Node.h"
#include "InOut.h"
#include "Integrator.h"
#include "InitialConditions.h"
#include "Analysis.h"


int main() {

//#pragma omp parallel
//	{
//		std::cout<<"test\n";
//	}
	//init stars
	int n_Stars = 3;
	double boxLength = 1; //[pc]
	double dt = 0.01;
	int nTimesteps = 2;
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

	std::vector<Star*> stars = {};

	//Init
	double totalMass = InitialConditions::initialMass(stars,n_Stars);
	InitialConditions::plummerSphere(stars, 1, totalMass);
	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;

	//Integrate
	Integrator rk4 = Integrator(dt);

	std::vector<double> time;
	std::vector<double> energy;
	std::vector<double> potE;
	std::vector<double> kinE;
	//Analysis::scaling(5, 5, euler);
	for (int i = 0; i < nTimesteps; i++) {

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::findCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr);
		for (Star* star : stars) {
			root.insert(star);
		}
		root.calculateMassDistribution();
		#pragma omp parallel for //1:10
		for (int i = 0; i < stars.size();++i){
			stars.at(i)->acceleration.reset();
			root.applyForce(stars.at(i));
		}
		rk4.RK4(stars,&root,dt);

		if (i % 10 == 0) {
			//InOut::writeWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::writeAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			double potentialEnergy = Analysis::potentialEnergy(stars);
			double kineticEnergy = Analysis::kineticEnergy(stars);
			std::cout << "stepnr: " << i << std::endl;
			std::cout<< "Kinetic Energy: " + std::to_string(kineticEnergy) << std::endl;
			std::cout << "Potential Energy: " + std::to_string(potentialEnergy) << std::endl;
			std::cout << "Total Energy: " + std::to_string(kineticEnergy+potentialEnergy) << std::endl << std::endl;
			time.push_back(i);
			energy.push_back(kineticEnergy + potentialEnergy);
			potE.push_back(potentialEnergy);
			kinE.push_back(kineticEnergy);
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(stars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	InOut::write(time, energy, "TotalEnergy.dat");
	InOut::write(time, kinE, "KinetikEnergy.dat");
	InOut::write(time, potE, "PotentialEnergy.dat");
	std::cout << "done" << std::endl;
	std::cin.get();
	return 0;
}