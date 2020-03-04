//#include <cstdlib>
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


int main() {

//#pragma omp parallel
//	{
//		std::cout<<"test\n";
//	}
	//init stars
	int n_Stars = 2;
	double boxLength = 1; //[pc]
	double dt = 1;
	int nTimesteps = 10;
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

	//std::random_device rd;  //Will be used to obtain a seed for the random number engine
	//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	//std::uniform_real_distribution<> dis(0.0, boxLength);
	std::vector<Star*> stars = {};

	//Init
	double totalMass = InitialConditions::InitialMass(stars,n_Stars);
	InitialConditions::PlummerSphere(stars, 1, totalMass);
	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;

	//Integrate
	Integrator euler = Integrator(dt);

	Analysis::scaling(10000, 10000, stars, euler); //WiP
	/*for (int i = 0; i < nTimesteps; i++) {

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::FindCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr);
		for (Star* star : stars) {
			root.Insert(star);
		}
		root.CalculateMassDistribution();
		#pragma omp parallel for //1:10
		for (int i = 0; i < stars.size();++i){
			stars.at(i)->acceleration.reset();
			root.ApplyForce(stars.at(i));
		}
		euler.Euler(stars,dt);

		if (i % 100 == 0) {
			InOut::WriteWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::WriteAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			//double potentialEnergy = Analysis::PotentialEnergy(stars);
			//double kineticEnergy = Analysis::KineticEnergy(stars);
			//std::cout<< "Kinetic Energy: " + std::to_string(kineticEnergy) << std::endl;
			//std::cout << "Potential Energy: " + std::to_string(potentialEnergy) << std::endl;
			//std::cout << "Total Energy: " + std::to_string(kineticEnergy+potentialEnergy) << std::endl << std::endl;
		}
	}
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(stars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;*/
	std::cout << "done" << std::endl;
	std::cin.get();
	return 0;
}