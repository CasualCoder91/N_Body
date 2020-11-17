#include "Analysis.h"

std::string Analysis::path = "src/Test/";

Analysis::Analysis(bool bEnergyDone, bool bVelocityDone, bool bVelocity2DDone){
	this->bEnergyDone = bEnergyDone;
	this->bVelocityDone = bVelocityDone;
	this->bVelocity2DDone = bVelocity2DDone;
}

double Analysis::potentialEnergy(std::vector<Star*>& stars){
	double potentialEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size()-1;++i) {
		for (int j = i+1; j < stars.size(); ++j) {
			potentialEnergy -= Constants::G * stars[i]->mass * stars[j]->mass / Vec3D::distance(&stars[i]->position, &stars[j]->position);
		}
	}
	this->potE.push_back(potentialEnergy);
	return potentialEnergy;
}

double Analysis::kineticEnergy(std::vector<Star*>& stars){
	double kineticEnergy = 0;
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars[i]->mass * (stars[i]->velocity * stars[i]->velocity);
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
		MWPotential potential = MWPotential();
		InitialConditions initialConditions = InitialConditions(&potential);
		int starID = 0;
		std::vector<Star*> stars = initialConditions.initStars(starID,Constants::nStars);
		double totalMass = initialConditions.initialMassSalpeter(stars, 0.08, 100);
		initialConditions.plummerSphere(stars, totalMass, Constants::boxLength, Constants::G);

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
				stars[i]->acceleration.reset(); // reset acceleration to 0,0,0
				root.applyForce(stars[i]);
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

bool Analysis::bTimestepAnalysis(){
	return this->bEnergyDone || this->bVelocityDone || this->bVelocity2DDone;
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

double Analysis::average(std::vector<Point>& points) {
	double average = 0;
	for (Point point : points) {
		average += sqrt(point.vx* point.vx + point.vy* point.vy);
	}
	return average / points.size();
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

double Analysis::dispersion(std::vector<Point>& points, double average){
	if (average == 0){ //average not passed as parameter
		average = Analysis::average(points);
	}

	double dispersion = 0;
	for (Point point : points) {
		dispersion += pow(sqrt(point.vx * point.vx + point.vy * point.vy) - average, 2);
	}
	return sqrt(dispersion / (points.size() - 1));
}

void Analysis::write(){
	if (time.size() != totE.size()) {
		std::cout << "time and energy vectors must have equal size! Aborting." << std::endl;
		return;
	}
	if (Constants::bEnergy) {
		InOut::write(time, totE, "TotalEnergy.dat");
		InOut::write(time, kinE, "KinetikEnergy.dat");
		InOut::write(time, potE, "PotentialEnergy.dat");
	}
}

void Analysis::cluster(std::vector<std::vector<Point>>& points){

	//for (int i = 0; i < points.size() - 1; ++i) { //loop through timesteps excluding last one
	for (Point& point0 : points[0]) { // loop through all points at timestep i
		double minDist = -1;
		Point futurePoint;
		for (Point& point1 : points[1]) { //compare to all points at timestep i+1
			double currentDist = point0.getDistance(point1);
			if (currentDist < minDist || minDist == -1) {
				minDist = currentDist;
				futurePoint = point1;
			}
		}
		//std::cout << point0.id << " " << futurePoint.id << std::endl;
		point0.vx = futurePoint.x-point0.x; //tecnically division by dt needed but dt is equal for all points
		point0.vy = futurePoint.y - point0.y;
	}
	//}

	//InOut::write(points[0], "src/Test/clusteringVelocity.dat");
	//Plot plot = Plot(path, path, true);
	//plot.plot("clusteringVelocity", {});

	VDBSCAN scanner = VDBSCAN(0.05, 0.0000001, 60);
	scanner.run(points[0]);

	int nFalsePositive = 0;
	int nFalseNegative = 0;
	int nTruePositive = 0;
	int nTrueNegative = 0;

	for (Point point : points[0]) {
		if (point.cluster > -1 && point.clusterStar)
			nTruePositive++;
		if (point.cluster > -1 && !point.clusterStar)
			nFalsePositive++;
		if (point.cluster == VDBSCAN::NOISE && point.clusterStar)
			nFalseNegative++;
		if (point.cluster == VDBSCAN::NOISE && !point.clusterStar)
			nTrueNegative++;
	}
	std::cout << "nFalsePositive: " << nFalsePositive << std::endl;
	std::cout << "nFalseNegative: " << nFalseNegative << std::endl;
	std::cout << "nTruePositive: " << nTruePositive << std::endl;
	std::cout << "nTrueNegative: " << nTrueNegative << std::endl;
}
