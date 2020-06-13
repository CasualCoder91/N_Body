#include "Analysis.h"

Analysis::Analysis(Parameters parameters){
	this->bEnergy = parameters.doEnergyAnalysis();
	this->bAverageVelocity = parameters.doAverageVelocity();
	this->bAverage2DVelocity = parameters.doAverage2DVelocity();
	this->G = parameters.getG();
}

double Analysis::potentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= this->G * stars.at(i)->mass * stars.at(j)->mass / Vec3D::distance(&stars.at(i)->position, &stars.at(j)->position);
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

void Analysis::scaling(int maxNStars, int nTimesteps, Integrator& integrator, Parameters* parameters){
	Vec3D tlf = Vec3D(), brb = Vec3D();

	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);

	std::vector<double> x, y;
	for (int n = 2; n <= maxNStars; ++n) {

		//init
		Potential potential = Potential(parameters);
		InitialConditions initialConditions = InitialConditions(&potential);
		int starID = 0;
		std::vector<Star*> stars = initialConditions.initStars(starID,parameters->getNStars());
		double totalMass = initialConditions.initialMassSalpeter(stars, 0.08, 100);
		initialConditions.plummerSphere(stars, totalMass,parameters->getBoxLength(),parameters->getG());

		startTime = std::chrono::steady_clock::now();
		for (int i = 0; i < nTimesteps; i++) {

			Node::findCorners(tlf, brb, stars);
			Node root = Node(tlf, brb, nullptr, parameters);
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

double Analysis::average(std::vector<Vec3D*>& vectors){
	double average = 0;
	for (Vec3D* vector : vectors) {
		average += vector->length();
	}
	return average/vectors.size();
}

double Analysis::average(std::vector<double>& values) {
	double average = 0;
	for (double value : values) {
		average += value;
	}
	return average / values.size();
}

double Analysis::dispersion(std::vector<Vec3D*>& vectors){
	size_t n = vectors.size();
	double average = Analysis::average(vectors);
	double dispersion = 0;

	for (Vec3D* vector : vectors) {
		dispersion += pow(vector->length()-average,2);
	}
	return sqrt(dispersion / (n-1));
}

double Analysis::dispersion(std::vector<double>& values) {
	size_t n = values.size();
	double average = Analysis::average(values);
	double dispersion = 0;

	for (double value : values) {
		dispersion += pow(value - average, 2);
	}
	return sqrt(dispersion / (n - 1));
}

void Analysis::write(){
	if (time.size() != totE.size()) {
		std::cout << "time and energy vectors must have equal size! Aborting." << std::endl;
		return;
	}
	if (bEnergy) {
		InOut::write(time, totE, "TotalEnergy.dat");
		InOut::write(time, kinE, "KinetikEnergy.dat");
		InOut::write(time, potE, "PotentialEnergy.dat");
	}
}
