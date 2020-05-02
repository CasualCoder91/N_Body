#include "InitialConditions.h"

double InitialConditions::kmInpc = 3.086e-13;

double InitialConditions::closestToZero(double a, double b) {
	double temp = 0;
	if ((a < 0 && b> 0) || (a > 0 && b < 0))
		temp = 0;
	if (a > 0 && b > 0) {
		if (a > b)
			return b;
		return a;
	}
	if (a < 0 && b < 0) {
		if (a > b)
			return a;
		return b;
	}
	return temp;
}

InitialConditions::InitialConditions(SimulationData* parameters):gen((std::random_device())()){
	this->G = parameters->getG();
	this->nStars = parameters->getNStars();
}

std::vector<Star*> InitialConditions::initStars(int firstID){
	std::vector<Star*> stars = {};
	for (int i = 0; i < this->nStars; i++) {
		Star* star = new Star(firstID+i);
		stars.push_back(star);
	}
	return stars;
}

std::vector<Star*> InitialConditions::initStars(int firstID, int nStars){
	std::vector<Star*> stars = {};
	for (int i = 0; i < nStars; i++) {
		Star* star = new Star(firstID + i);
		stars.push_back(star);
	}
	return stars;
}

double InitialConditions::initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha){
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0);//avoid close to singularity
	double totalMass = 0.0;
	double temp = alpha + 1.0;
	double factor = (pow(maxMass / minMass, temp) - 1.0);
	for (Star* star : stars) {
		star->mass = minMass * pow(1.0 + factor * dis(gen), 1.0 / temp);
		totalMass += star->mass;
	}
	return totalMass;
}

std::vector<Star*> InitialConditions::initDiskStars(int firstID, Vec3D tlf, Vec3D brf, double depth, Potential* potential, double gridResolution){
	std::vector<Star*> stars;
	if (depth < 0) {
		throw "initDiskStars: Depth must be >= 0";
	}
	for (double z = tlf.z; z < tlf.z + depth; z += gridResolution) {
		std::cout << "z = " << z << " max z = " << tlf.z + depth << std::endl;
		for (double x = tlf.x; x < brf.x; x += gridResolution) {
			for (double y = brf.y; y < tlf.y; y += gridResolution) {
				Vec3D position = Vec3D(x, y, z);
				Vec3D volumeElement = Vec3D(gridResolution, gridResolution, gridResolution);
				double massInCell = potential->massDisk(position, volumeElement);
				std::vector<Star*> starsInCell = massDisk(massInCell); //stars with mass
				sampleDiskPositions(starsInCell, position, volumeElement);
				stars.insert(stars.end(), starsInCell.begin(), starsInCell.end());
			}
		}
	}

	return stars;
}

double InitialConditions::sampleDiskPositions(std::vector<Star*> stars,Vec3D position, Vec3D volumeElement) {
	double smallestx = closestToZero(position.x, position.x+volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);

	double k = Potential::densityDisk(smallestx, smallesty, smallestz);
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(position.x, position.x + volumeElement.x);
	std::uniform_real_distribution<> disy(position.y, position.y + volumeElement.y);
	std::uniform_real_distribution<> disz(position.z, position.z + volumeElement.z);
	std::uniform_real_distribution<> disaccept(0, k);
	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);

			double accept = disaccept(gen);
			double temp = Potential::densityDisk(x, y, z);

			if (accept < temp) {
				if(sqrt(pow(x,2) + pow(y,2)) < 2000)
					continue;
				star->position = Vec3D(x, y, z);
				break;
			}
		}
	}
	return 0; //todo: return average velocity maybe?
}

void InitialConditions::sampleDiskVelocity(Vec3D& velocity, Vec3D& position){

	double R = gsl_hypot(position.x,position.y);

	std::normal_distribution<> zVelocityDistribution{ 0,Potential::verticalVelocityDispersion(R) };
	double vz = zVelocityDistribution(gen);

	double rDispersion = Potential::radialVelocityDispersionDisk(R,position.z);
	double vR = 0;
	if (rDispersion > 0){
		std::normal_distribution<> radialVelocityDistribution{ 0,rDispersion };
		vR = radialVelocityDistribution(gen);
	}
	double aDispersion = Potential::azimuthalVelocityDispersion(R, position.z);
	double va = 0;
	if (aDispersion > 0) {
		std::normal_distribution<> azimuthalVelocityDistribution{ Potential::azimuthalStreamingVelocity(position),aDispersion };
		va = azimuthalVelocityDistribution(gen);
	}
	else {
		va = Potential::azimuthalStreamingVelocity(position);
	}

	double theta = atan2(position.y , position.x);

	velocity += Vec3D(vR * cos(theta) + va * cos(theta+M_PI_2), vR * sin(theta) + va * sin(theta + M_PI_2), vz);
}

double InitialConditions::sampleDiskVelocities(std::vector<Star*> stars){
	//std::random_device rd{};
	//std::mt19937 gen{ rd() };
	for (Star* star : stars) {
		sampleDiskVelocity(star->velocity, star->position);
	}
	return 0.0;
}

double InitialConditions::sampleBulgePositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement){
	double smallestx = closestToZero(position.x, position.x + volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);

	double k = Potential::densityDisk(smallestx, smallesty, smallestz);
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(position.x, position.x + volumeElement.x);
	std::uniform_real_distribution<> disy(position.y, position.y + volumeElement.y);
	std::uniform_real_distribution<> disz(position.z, position.z + volumeElement.z);
	std::uniform_real_distribution<> disaccept(0, k);
	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);

			double accept = disaccept(gen);
			double temp = Potential::densityBulge(x, y, z);

			if (accept < temp) {
				star->position = Vec3D(x, y, z);
				break;
			}
		}
	}
	return 0; //todo: return average velocity maybe?
}

//void InitialConditions::sampleBulgeVelocity(Vec3D& velocity, Vec3D& position){
//	//std::random_device rd;  //Will be used to obtain a seed for the random number engine
//	//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//	//double escapeVelocity = sqrt(-2 * Potential::potentialEnergy(position));
//	//std::uniform_real_distribution<> distributionVelocity(0.0, escapeVelocity); // Use escape velocity instead!
//	//double maximumEnergy = Potential::particleEnergy(position, Vec3D(0, 0, 0));
//	//double maxDF = Potential::distributionFunctionBulge(-maximumEnergy);//(gsl_pow_2(Potential::characteristicVelocityBulge)-1); // -1 to avoid division by zero
//	//std::uniform_real_distribution<> distributionDF(0.0, maxDF);
//	//while (true) {
//	//	double mVelocity = distributionVelocity(gen);
//	//	double randomDF = distributionDF(gen);
//	//	double energy = Potential::particleEnergy(position, mVelocity);
//	//	double attemptDF = Potential::distributionFunctionBulge(-energy);
//	//	if (randomDF < attemptDF) {
//	//		position += Vec3D::randomVector(mVelocity);
//	//		return;
//	//	}
//	//}
//	double R = gsl_hypot(position.x, position.y);
//	//std::cout << position.print() << std::endl;
//	std::normal_distribution<> velocityDistribution{ 0,Potential::radialVelocityDispersionBulge(R, position.z) };
//	double vz = velocityDistribution(gen);
//	double vR = velocityDistribution(gen);
//	//std::normal_distribution<> velocityDistribution2{ 200,2*Potential::radialVelocityDispersionBulge(R, position.z) };
//	double va = velocityDistribution(gen);
//
//	double theta = atan2(position.y, position.x);
//
//	velocity += Vec3D(vR * cos(theta) + va * cos(theta + M_PI_2), vR * sin(theta) + va * sin(theta + M_PI_2), vz);
//}

std::vector<Star*> InitialConditions::massDisk(double totalMass){
	std::vector<Star*> stars;
	if (totalMass <= 0)
		return stars;
	double pickedTotalMass = 0;
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> disM(0.08, 63.1); // mass sample
	static std::uniform_real_distribution<> disAccept(0, 0.86); //upper limit
	static double factor1 = 0.158;
	static double chabrierMass = 0.079;
	static double chabrierSigma = 0.69;
	static double factor2 = 4.4e-2;
	static double exponent2 = 4.37;
	static double factor3 = 1.5e-2;
	static double exponent3 = 3.53;
	static double factor4 = 2.5e-4;
	static double exponent4 = 2.11;
	double temp = 0;
	double ln10 = log(10);
	while (pickedTotalMass < totalMass) {
		while (true) {
			double m = disM(gen);
			double logM = log10(m);
			if (logM < 0) { // m < 1
				temp = factor1/(m* ln10)* exp(-pow(logM - log10(chabrierMass), 2) / (2.0 * pow(chabrierSigma, 2)));
				if (disAccept(gen) < temp ) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					break;
				}
			}
			else if (logM < 0.54) {
				temp = factor2 / (m * ln10) * pow(m, -exponent2);
				if (disAccept(gen) < temp) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					break;
				}
			}
			else if (logM < 1.26){
				temp = factor3 / (m * ln10) * pow(m, -exponent3);
				if (disAccept(gen)< temp) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					std::cout << "pickedTotalMass: " << pickedTotalMass << std::endl;
					break;
				}
			}
			else{
				temp = factor4 / (m * ln10) * pow(m, -exponent4);
				if (disAccept(gen) < temp) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					break;
				}
			}
		}
	}
	std::cout << "Mean mass: " << pickedTotalMass / stars.size() << std::endl;
	std::cout << "Proposed mass of disk: " << totalMass << " Sampled mass: " << pickedTotalMass << std::endl;
	return stars;
}


std::vector<Star*> InitialConditions::initialMassBulge(double totalMass){
	std::vector<Star*> stars;
	double pickedTotalMass = 0;
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> disM(0.08, 1); // mass sample
	static std::uniform_real_distribution<> disAccept(0, 0.00095); //upper limit
	static double factor1 = 3.6e-4;
	static double chabrierMass = 0.22;
	static double chabrierSigma = 0.33;
	static double factor2 = 7.1e-5;
	static double exponent2 = 1.3;
	double temp = 0;
	double ln10 = log(10);
	while (pickedTotalMass < totalMass) {
		while (true) {
			double m = disM(gen);
			double logM = log10(m);
			if (m<0.7) { // m < 1
				temp = factor1 / (m * ln10) * exp(-pow(logM - log10(chabrierMass), 2) / (2.0 * pow(chabrierSigma, 2)));
				if (disAccept(gen) < temp) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					break;
				}
			}
			else{
				temp = factor2 / (m * ln10) * pow(m, -exponent2);
				if (disAccept(gen) < temp) {
					stars.push_back(new Star(0, m)); // todo: which id?!
					pickedTotalMass += m;
					break;
				}
			}
		}
	}
	std::cout << "Mean mass: " << pickedTotalMass / stars.size() << std::endl;
	std::cout << "Proposed mass of disk: " << totalMass << " Sampled mass: " << pickedTotalMass << std::endl;
	return stars;
}

void InitialConditions::plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass){
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 0.99);//avoid close to singularity
	for (Star* star : stars) {
		double distance = structuralLength / sqrt(pow(dis(gen), -2. / 3.) - 1);
		star->position = Vec3D::randomVector(distance);
		plummerVelocity(star, structuralLength, distance, totalMass);
	}
}

void InitialConditions::offsetCluster(std::vector<Star*>& stars, Vec3D& offset){
	for (Star* star : stars) {
		star->position += offset;
	}
}

void InitialConditions::setNStars(int N){
	this->nStars = N;
}

double InitialConditions::plummerEscapeVelocity(double distance, double structuralLength, double totalMass){
	//return sqrt(2.) * pow(distance * distance + structuralLength, -0.25);
	//https://github.com/bacook17/behalf/blob/master/behalf/initialConditions.py
	return sqrt(2. * this->G * totalMass / structuralLength) *pow(1.+distance*distance/(structuralLength* structuralLength),-0.25);
}

void InitialConditions::plummerVelocity(Star* star, double structuralLength, double distance, double totalMass){
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0);
	//Rejection Technique
	std::uniform_real_distribution<> dis_g(0.0, 0.1);
	double q = dis(gen); // random value in range [0,1]
	double g = dis_g(gen); // random value in range [0,0.1]
	while (g > q* q* pow( (1. - q * q), 3.5)) {
		q = dis(gen);
		g = dis_g(gen);
	}
	double velocity = q * plummerEscapeVelocity(distance, structuralLength, totalMass)/1.01;
	star->velocity = Vec3D::randomAngles(velocity);
}
