#include "Analysis.h"

double Analysis::PotentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	for (unsigned int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= Parameters::G * stars.at(i)->mass * stars.at(j)->mass / Vec3D::distance(&stars.at(i)->position, &stars.at(j)->position);
		}
	}
	return potentialEnergy;
}

double Analysis::KineticEnergy(std::vector<Star*>& stars){
	double kineticEnergy = 0;
	for (unsigned int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars.at(i)->mass * (stars.at(i)->velocity * stars.at(i)->velocity);
	}
	kineticEnergy *= 0.5;
	return kineticEnergy;
}

void Analysis::scaling(int maxN, int nTimesteps, std::vector<Star*>& stars, Integrator& integrator){
	Vec3D tlf = Vec3D(), brb = Vec3D();
	Node root = Node(tlf, brb, nullptr);
	std::vector<double> x, y;
	for (int n = 0; n < maxN; n++) {
		for (int i = 0; i < nTimesteps; i++) {
			std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

			Node::FindCorners(tlf, brb, stars);
			root = Node(tlf, brb, nullptr);
			for (Star* star : stars) {
				root.Insert(star);
			}
			root.CalculateMassDistribution();
			#pragma omp parallel for //1:10
			for (int i = 0; i < stars.size(); ++i) {
				stars.at(i)->acceleration.reset(); // reset acceleration to 0,0,0
				root.ApplyForce(stars.at(i));
			}
			integrator.Euler(stars);

			std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
			x.push_back(n);
			y.push_back(time_span.count());

			//std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
			//std::cout << "done" << std::endl;
			//std::cin.get();
		}
	}
	InOut::write(x,y,"NlogNtest.dat");
}
