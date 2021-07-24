#include "Analysis.h"

std::string Analysis::path = "src/Test/";

Analysis::Analysis(int id, Database* database)
{
	this->id = id;
	this->database = database;
}


double Analysis::potentialEnergy(const std::vector<Star>& stars){
	double potentialEnergy = 0;
    #pragma omp parallel for reduction(+:potentialEnergy)
	for (int i = 0; i < stars.size()-1;++i) {
		//for (int j = i+1; j < stars.size(); ++j) {
			//potentialEnergy -= Constants::G * stars[i].mass * stars[j].mass / Vec3D::distance(stars[i].position, stars[j].position);
		//}
		double temp = stars[i].mass* MWPotential::potentialEnergy(stars[i].position);
		potentialEnergy = potentialEnergy + temp;
	}
	potentialEnergy *= 0.5; //sure?
	this->potE.push_back(potentialEnergy);
	return potentialEnergy;
}

double Analysis::kineticEnergy(const std::vector<Star>& stars){
	double kineticEnergy = 0;
    #pragma omp parallel for reduction(+:kineticEnergy)
	for (int i = 0; i < stars.size(); ++i) {
		kineticEnergy += stars[i].mass *(stars[i].velocity.x* stars[i].velocity.x + stars[i].velocity.y* stars[i].velocity.y + stars[i].velocity.z* stars[i].velocity.z);
	}
	kineticEnergy *= 0.5;
	this->kinE.push_back(kineticEnergy);
	return kineticEnergy;
}

void Analysis::energy()
{
	std::vector<int> timeSteps = database->selectTimesteps(id);
	for (int timeStep : timeSteps) {
		std::vector<Star> stars = database->selectStars(id, timeStep);
		double kinE = kineticEnergy(stars);
		double potE = potentialEnergy(stars);
		database->insertAnalysisdtEnergy(id, timeStep, kinE, potE);
	}
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

double Analysis::average(std::vector<Vec3D>& vectors){
	double average = 0;
	for (Vec3D vector : vectors) {
		average += vector.length();
	}
	return average/vectors.size();
}

double Analysis::average(std::vector<Vec2D>& vectors) {
	double average = 0;
	for (Vec2D vector : vectors) {
		average += vector.length();
	}
	return average / vectors.size();
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
		average += sqrt(point.velocity[0]* point.velocity[0] + point.velocity[1]* point.velocity[1]);
	}
	return average / points.size();
}

double Analysis::dispersion(std::vector<Vec3D>& vectors, double average){
	size_t n = vectors.size();
	if (average == -1) { //average not passed as parameter
		average = Analysis::average(vectors);
	}

	double dispersion = 0;
	for (Vec3D vector : vectors) {
		dispersion += pow(vector.length()-average,2);
	}
	return sqrt(dispersion / (n-1));
}

double Analysis::dispersion(std::vector<Vec2D>& vectors, double average) {
	size_t n = vectors.size();
	if (average == -1) { //average not passed as parameter
		average = Analysis::average(vectors);
	}

	double dispersion = 0;
	for (Vec2D vector : vectors) {
		dispersion += pow(vector.length() - average, 2);
	}
	return sqrt(dispersion / (n - 1));
}

double Analysis::dispersion(std::vector<double>& values, double average) {
	size_t n = values.size();
	if (average == -1) { //average not passed as parameter
		average = Analysis::average(values);
	}
	double dispersion = 0;

	for (double value : values) {
		dispersion += pow(value - average, 2);
	}
	return sqrt(dispersion / (n - 1));
}

double Analysis::dispersion(std::vector<Point>& points, double average){
	if (average == -1){ //average not passed as parameter
		average = Analysis::average(points);
	}

	double dispersion = 0;
	for (Point point : points) {
		dispersion += pow(sqrt(point.velocity[0] * point.velocity[0] + point.velocity[1] * point.velocity[1]) - average, 2);
	}
	return sqrt(dispersion / (points.size() - 1));
}

void Analysis::write(){
	if (time.size() != totE.size()) {
		std::cout << "time and energy vectors must have equal size! Aborting." << std::endl;
		return;
	}
	if (totE.size()>0) {
		InOut::write(time, totE, "TotalEnergy.dat");
		InOut::write(time, kinE, "KinetikEnergy.dat");
		InOut::write(time, potE, "PotentialEnergy.dat");
	}
}

void Analysis::generateHTPVelocity(bool observed)
{
	std::vector<std::vector<Point>>& points = database->selectPoints(id, 0, 2,-1,observed);

	int errorCounter = 0;

	for (Point& point0 : points[0]) { // loop through all points at timestep i
		double minDist = -1;
		Point futurePoint;
		for (Point& point1 : points[1]) { //compare to all points at timestep i+1
			double currentDist = point0.getDistance(point1);
			if ((currentDist < minDist || minDist == -1) && abs(1 - point1.magnitude / point0.magnitude) < Constants::epsMagnitude) {
				minDist = currentDist;
				futurePoint = point1;
			}
		}
		if (point0.id != futurePoint.id) {
			errorCounter++;
		}
		point0.velocity[0] = futurePoint.x - point0.x; //tecnically division by dt needed but dt is equal for all points
		point0.velocity[1] = futurePoint.y - point0.y;
	}
	std::cout << "#Wrong stars picked for velocity estimation: " << errorCounter << std::endl;
	database->updatePoints(points[0]);
}


void Analysis::cluster(std::vector<Point>& points) {

	mlpack::dbscan::DBSCAN<> dbscan(7.46105964516358e-06*0.3, 200);

	arma::mat matPoints = arma::mat(2, points.size());
	for (size_t i = 0; i < points.size();i++) {
		arma::vec column = arma::vec(2);
		column.at(0) = points[i].velocity[0];
		column.at(1) = points[i].velocity[1];
		matPoints.col(i) = column;
	}
	arma::Row<size_t> assignments;
	const size_t nClusters = dbscan.Cluster(matPoints, assignments);

	printf("%zu cluster(s) detected\n", nClusters);

	int nFalsePositive = 0;
	int nFalseNegative = 0;
	int nTruePositive = 0;
	int nTrueNegative = 0;

	for (size_t i = 0; i < points.size(); i++) {
		if (assignments[i] == SIZE_MAX)
			points[i].cluster = -1;
		else
			points[i].cluster = assignments[i];
		if(points[i].cluster == 0 && points[i].clusterStar)
			nTruePositive++;
		else if(points[i].cluster != 0 && points[i].clusterStar)
			nFalseNegative++;
		else if (points[i].cluster == -1 && !points[i].clusterStar)
			nTrueNegative++;
		else if (points[i].cluster != -1 && !points[i].clusterStar)
			nFalsePositive++;
	}

	std::cout << "nFalsePositive: " << nFalsePositive << std::endl;
	std::cout << "nFalseNegative: " << nFalseNegative << std::endl;
	std::cout << "nTruePositive: " << nTruePositive << std::endl;
	std::cout << "nTrueNegative: " << nTrueNegative << std::endl;
	//points.push_back(resultAt0);
	database->updatePoints(points);
}
