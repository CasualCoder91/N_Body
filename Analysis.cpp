#include "Analysis.h"

double Analysis::potentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= Parameters::G * stars.at(i)->mass * stars.at(j)->mass / Vec3D::distance(&stars.at(i)->position, &stars.at(j)->position);
		}
	}
	return potentialEnergy;
}

double Analysis::kineticEnergy(std::vector<Star*>& stars){
	double kineticEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars.at(i)->mass * (stars.at(i)->velocity * stars.at(i)->velocity);
	}
	kineticEnergy *= 0.5;
	return kineticEnergy;
}

void Analysis::scaling(int maxNStars, int nTimesteps, Integrator& integrator){
	Vec3D tlf = Vec3D(), brb = Vec3D();

	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);

	std::vector<double> x, y;
	for (int n = 2; n <= maxNStars; ++n) {

		//init
		std::vector<Star*> stars = {};
		double totalMass = InitialConditions::initialMass(stars, n);
		InitialConditions::plummerSphere(stars, 1, totalMass);

		startTime = std::chrono::steady_clock::now();
		for (int i = 0; i < nTimesteps; i++) {

			Node::findCorners(tlf, brb, stars);
			Node root = Node(tlf, brb, nullptr);
			for (Star* star : stars) {
				root.insert(star);
			}
			root.calculateMassDistribution();
			#pragma omp parallel for //1:10
			for (int i = 0; i < stars.size(); ++i) {
				stars.at(i)->acceleration.reset(); // reset acceleration to 0,0,0
				root.applyForce(stars.at(i));
			}
			integrator.euler(stars);

			//std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
			//std::cout << "done" << std::endl;
			//std::cin.get();
		}
		endTime = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
		x.push_back(n);
		y.push_back(time_span.count());
		std::cout << "n: " << n << " N: " << maxNStars << std::endl;
	}
	InOut::write(x,y,"Output/NlogNtest.dat");
}
