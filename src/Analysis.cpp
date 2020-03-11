#include "Analysis.h"

Analysis::Analysis(){
	Parameters::config.Load(filePath);
	if (!config.Get("bEnergy", bEnergy)) {
		bEnergy = false;
		std::cout << "bEnergy missing in " + filePath + "add bEnergy = true to activate energy analysis."<< std::endl;
	}
	if (!config.Get("bAverageVelocity", bAverageVelocity)) {
		bAverageVelocity = false;
		std::cout << "bAverageVelocity missing in " + filePath + "add bAverageVelocity = true to activate 3d velocity analysis." << std::endl;
	}
	if (!config.Get("bAverage2DVelocity", bAverage2DVelocity)) {
		bAverageVelocity = false;
		std::cout << "bAverage2DVelocity missing in " + filePath + "add bAverage2DVelocity = true to activate 2d velocity analysis." << std::endl;
	}
}

double Analysis::potentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= Parameters::G * stars.at(i)->mass * stars.at(j)->mass / Vec3D::distance(&stars.at(i)->position, &stars.at(j)->position);
		}
	}
	this->potE.push_back(potentialEnergy);
	return potentialEnergy;
}

double Analysis::kineticEnergy(std::vector<Star*>& stars){
	double kineticEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars.at(i)->mass * (stars.at(i)->velocity * stars.at(i)->velocity);
	}
	kineticEnergy *= 0.5;
	this->kinE.push_back(kineticEnergy);
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
		InitialConditions initialConditions = InitialConditions();
		std::vector<Star*> stars = initialConditions.initStars(0,n);
		double totalMass = initialConditions.initialMass(stars);
		initialConditions.plummerSphere(stars, 1, totalMass);

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
			integrator.euler(stars,&root);
		}
		endTime = std::chrono::steady_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
		x.push_back(n);
		y.push_back(time_span.count());
		std::cout << "n: " << n << " N: " << maxNStars << std::endl;
	}
	InOut::write(x,y,"NlogNtest.dat");
}

bool Analysis::getbEnergy(){
	return this->bEnergy;
}

bool Analysis::getbAverageVelocity(){
	return this->bAverageVelocity;
}

bool Analysis::getbAverage2DVelocity(){
	return this->bAverage2DVelocity;
}

void Analysis::write(){
	if (bEnergy) {
		InOut::write(time, totE, "TotalEnergy.dat");
		InOut::write(time, kinE, "KinetikEnergy.dat");
		InOut::write(time, potE, "PotentialEnergy.dat");
	}
}
