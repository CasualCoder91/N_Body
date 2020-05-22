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

double InitialConditions::farthermostFromZero(double a, double b){
	double temp = 0;
	if (abs(a) > abs(b))
		return a;
	return b;
}

InitialConditions::InitialConditions(Potential* potential):gen((std::random_device())()){
	this->potential = potential;
}

std::vector<Star*> InitialConditions::initStars(int& firstID,int nStars){
	std::vector<Star*> stars = {};
	for (; firstID < nStars; firstID++) {
		Star* star = &Star(firstID);
		stars.push_back(star);
	}
	return stars;
}

std::vector<Star*> InitialConditions::initFieldStars(int& starID, Vec3D focus, Vec3D viewPoint, double distance, double dx, double angle){
	std::cout << "Initializing field stars" << std::endl;
	angle = angle * 0.0174533; //convert degrees in rad
	Vec3D direction = (focus - viewPoint).normalize();
	int nSteps = (distance+dx) / dx;
	ProgressBar progressBar = ProgressBar(0, nSteps, true);
	std::vector<Star*> fieldStars;
	//#pragma omp parallel for
	for (int step = 1; step <= nSteps; ++step) {//steps along direction (line of sight)
		double r = 0.5*step*dx * tan(angle); //distance from line of sight at current step
		double aBoid = 2 * r;
		if (aBoid < 1) //if cubes are smaller 1pc^3 density is aproximated 0
			continue;
		Vec3D rVec = sqrt(2) * r * Vec3D::crossProduct(&Vec3D(-1, 1, -1), &direction).normalize();
		rVec.x = -abs(rVec.x);
		rVec.y = -abs(rVec.y);
		rVec.z = -abs(rVec.z);
		Vec3D corner = viewPoint + direction * ((double)step - 1)*dx + rVec;
		//std::cout << corner.print() << std::endl;
		Vec3D volumeElement = direction*dx-2*rVec;
		double diskMass = potential->massDisk(corner, volumeElement);
		std::vector<Star*> diskStars = massDisk(diskMass,starID);
		if (diskStars.size() > 0) {
			sampleDiskPositions(diskStars, corner, volumeElement); //test
			sampleDiskVelocities(diskStars);
			fieldStars.insert(std::end(fieldStars), std::begin(diskStars), std::end(diskStars));
		}
		double bulgeMass = potential->massBulge(corner, volumeElement);
		std::vector<Star*> bulgeStars = initialMassBulge(bulgeMass,starID);
		if (bulgeStars.size() > 0) {
			sampleBulgePositions(bulgeStars, corner, volumeElement); //test
			//sampleBulgeVelocities(bulgeStars);
			fieldStars.insert(std::end(fieldStars), std::begin(bulgeStars), std::end(bulgeStars));
		}
		progressBar.Update(step);
		progressBar.Print();
	}
	return fieldStars;
}

double InitialConditions::bulgeStarMass(Vec3D focus, Vec3D viewPoint, double distance, double dx, double angle){
	std::cout << "Calculating Bulge Star Mass" << std::endl;
	angle = angle * 0.0174533; //convert degrees in rad
	Vec3D direction = (focus - viewPoint).normalize();
	int nSteps = (distance + dx) / dx;
	ProgressBar progressBar = ProgressBar(0, nSteps, true);
	double totalMass = 0;
	//#pragma omp parallel for
	for (int step = 1; step <= nSteps; ++step) {//steps along direction (line of sight)
		double r = 0.5 * step * dx * tan(angle); //distance from line of sight at current step
		double aBoid = 2 * r;
		if (aBoid < 1) //if cubes are smaller 1pc^3 density is aproximated 0
			continue;
		Vec3D rVec = sqrt(2) * r * Vec3D::crossProduct(&Vec3D(-1, 1, -1), &direction).normalize();
		rVec.x = -abs(rVec.x);
		rVec.y = -abs(rVec.y);
		rVec.z = -abs(rVec.z);
		Vec3D corner = viewPoint + direction * ((double)step - 1) * dx + rVec;
		//std::cout << corner.print() << std::endl;
		Vec3D volumeElement = direction * dx - 2 * rVec;
		totalMass += potential->massBulge(corner, volumeElement);
		progressBar.Update(step);
		progressBar.Print();
	}
	return totalMass;
}

double InitialConditions::diskStarMass(Vec3D focus, Vec3D viewPoint, double distance, double dx, double angle){
	std::cout << "Calculating Disk Star Mass" << std::endl;
	angle = angle * 0.0174533; //convert degrees in rad
	Vec3D direction = (focus - viewPoint).normalize();
	int nSteps = (distance + dx) / dx;
	ProgressBar progressBar = ProgressBar(0, nSteps, true);
	double totalMass = 0;
	//#pragma omp parallel for
	for (int step = 1; step <= nSteps; ++step) {//steps along direction (line of sight)
		double r = 0.5 * step * dx * tan(angle); //distance from line of sight at current step
		double aBoid = 2 * r;
		if (aBoid < 1) //if cubes are smaller 1pc^3 density is aproximated 0
			continue;
		Vec3D rVec = sqrt(2) * r * Vec3D::crossProduct(&Vec3D(-1, 1, -1), &direction).normalize();
		rVec.x = -abs(rVec.x);
		rVec.y = -abs(rVec.y);
		rVec.z = -abs(rVec.z);
		Vec3D corner = viewPoint + direction * ((double)step - 1) * dx + rVec;
		if (corner.length() > potential->aBulge)
			continue;
		//std::cout << corner.print() << std::endl;
		Vec3D volumeElement = direction * dx - 2 * rVec;
		totalMass += potential->massDisk(corner, volumeElement);
		progressBar.Update(step);
		progressBar.Print();
	}
	return totalMass;
}

double InitialConditions::initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha){
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

std::vector<Star*> InitialConditions::initDiskStars(int& starID, Vec3D tlf, Vec3D brf, double depth, double gridResolution){
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
				std::vector<Star*> starsInCell = massDisk(massInCell, starID); //stars with mass
				sampleDiskPositions(starsInCell, position, volumeElement);
				stars.insert(stars.end(), starsInCell.begin(), starsInCell.end());
			}
		}
	}

	return stars;
}

double InitialConditions::sampleDiskPositions(std::vector<Star*> stars,Vec3D position, Vec3D volumeElement) {
	//get Limits for "Accept" distribution. Depends on Volume. Largest density is found at positions closest to 0
	double smallestx = closestToZero(position.x, position.x+volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);
	double largestx = farthermostFromZero(position.x, position.x + volumeElement.x);
	double largesty = farthermostFromZero(position.y, position.y + volumeElement.y);
	double largestz = farthermostFromZero(position.z, position.z + volumeElement.z);
	double acceptUpperLimit = potential->densityDisk(smallestx, smallesty, smallestz);
	double acceptLowerLimit = potential->densityDisk(largestx, largesty, largestz);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);//lower limit is new -> speedup

	double x1_low = std::min(position.x, position.x + volumeElement.x);
	double x1_high = std::max(position.x, position.x + volumeElement.x);
	double x2_low = std::min(position.y, position.y + volumeElement.y);
	double x2_high = std::max(position.y, position.y + volumeElement.y);
	double x3_low = std::min(position.z, position.z + volumeElement.z);
	double x3_high = std::max(position.z, position.z + volumeElement.z);
	std::uniform_real_distribution<> disx(x1_low, x1_high);
	std::uniform_real_distribution<> disy(x2_low, x2_high);
	std::uniform_real_distribution<> disz(x3_low, x3_high);

	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);

			double accept = disaccept(gen);
			double temp = potential->densityDisk(x, y, z);

			if (accept < temp) {
				star->position = Vec3D(x, y, z);
				break;
			}
		}
	}
	return 0; //todo: return average velocity maybe?
}

void InitialConditions::sampleDiskVelocity(Vec3D& velocity, Vec3D& position){

	double R = gsl_hypot(position.x,position.y);

	std::normal_distribution<> zVelocityDistribution{ 0,potential->verticalVelocityDispersion(R) };
	double vz = zVelocityDistribution(gen);

	double rDispersion = potential->radialVelocityDispersionDisk(R,position.z);
	double vR = 0;
	if (rDispersion > 0){
		std::normal_distribution<> radialVelocityDistribution{ 0,rDispersion };
		vR = radialVelocityDistribution(gen);
	}
	double aDispersion = potential->azimuthalVelocityDispersion(R, position.z);
	double va = 0;
	if (aDispersion > 0) {
		std::normal_distribution<> azimuthalVelocityDistribution{ potential->azimuthalStreamingVelocity(position),aDispersion };
		va = azimuthalVelocityDistribution(gen);
	}
	else {
		va = potential->azimuthalStreamingVelocity(position);
	}

	double theta = atan2(position.y , position.x);

	velocity += Vec3D(vR * cos(theta) + va * cos(theta+M_PI_2), vR * sin(theta) + va * sin(theta + M_PI_2), vz);
}

double InitialConditions::sampleDiskVelocities(std::vector<Star*> stars){
	for (Star* star : stars) {
		sampleDiskVelocity(star->velocity, star->position);
	}
	return 0.0;
}

void InitialConditions::sampleBulgePositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement){
	double smallestx = closestToZero(position.x, position.x + volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);
	double largestx = farthermostFromZero(position.x, position.x + volumeElement.x);
	double largesty = farthermostFromZero(position.y, position.y + volumeElement.y);
	double largestz = farthermostFromZero(position.z, position.z + volumeElement.z);
	double acceptUpperLimit = potential->densityBulge(smallestx, smallesty, smallestz);
	double acceptLowerLimit = potential->densityBulge(largestx, largesty, largestz);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);//lower limit is new -> speedup


	double x1_low = std::min(position.x, position.x + volumeElement.x);
	double x1_high = std::max(position.x, position.x + volumeElement.x);
	double x2_low = std::min(position.y, position.y + volumeElement.y);
	double x2_high = std::max(position.y, position.y + volumeElement.y);
	double x3_low = std::min(position.z, position.z + volumeElement.z);
	double x3_high = std::max(position.z, position.z + volumeElement.z);
	//double rMin = gsl_hypot3(std::min(position.x, position.x + volumeElement.x), std::min(position.y, position.y + volumeElement.y), std::min(position.z, position.z + volumeElement.z));
	//double rMax = gsl_hypot3(std::max(position.x, position.x + volumeElement.x), std::max(position.y, position.y + volumeElement.y), std::max(position.z, position.z + volumeElement.z));
	//std::uniform_real_distribution<> disr(rMin, rMax);
	std::uniform_real_distribution<> disx(x1_low, x1_high);
	std::uniform_real_distribution<> disy(x2_low, x2_high);
	std::uniform_real_distribution<> disz(x3_low, x3_high);

	//const double lamda = 1 / (Potential::aBulge * 2);
	//std::exponential_distribution<double> disr(lamda);
	//std::exponential_distribution<double> disaccept(lamda);
	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);
			double accept = disaccept(gen);
			double temp = potential->densityBulge(x, y, z);

			if (accept < temp) {
				star->position = Vec3D(x, y, z);
				break;
			}
		}
	}
	return; //todo: return average velocity maybe?
}

void InitialConditions::sampleBulgeVelocity(Vec3D& velocity, Vec3D& position){
	double delta = potential->velocityDistributionBulgeTableValue(position.length());
	double vCirc = potential->circularVelocity(&position);
	std::normal_distribution<> velocityDistribution{ 0,delta };
	double vRand = velocityDistribution(gen);

	double theta = atan2(position.y, position.x);
	double phi = asin(position.z / position.length());
	velocity = Vec3D(vCirc * cos(theta + M_PI_2) + cos(phi)*cos(theta)*vRand, vCirc * sin(theta + M_PI_2) + cos(phi) * sin(theta) * vRand, sin(phi)*vRand);

}

void InitialConditions::sampleBulgeVelocities(std::vector<Star*> stars){
	for (Star* star : stars) {
		sampleBulgeVelocity(star->velocity, star->position);
	}
}

std::vector<Star*> InitialConditions::massDisk(double totalMass, int& starID){
	std::vector<Star*> stars;
	if (totalMass <= 0)
		return stars;
	Star* proposedStar;
	double pickedTotalMass = 0;
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
					proposedStar = new Star(0, m);
					break;
				}
			}
			else if (logM < 0.54) {
				temp = factor2 / (m * ln10) * pow(m, -exponent2);
				if (disAccept(gen) < temp) {
					proposedStar = new Star(0, m);
					break;
				}
			}
			else if (logM < 1.26){
				temp = factor3 / (m * ln10) * pow(m, -exponent3);
				if (disAccept(gen)< temp) {
					proposedStar = new Star(0, m);
					break;
				}
			}
			else{
				temp = factor4 / (m * ln10) * pow(m, -exponent4);
				if (disAccept(gen) < temp) {
					proposedStar = new Star(0, m);
					break;
				}
			}
		}
		if (pickedTotalMass + proposedStar->mass / 2 < totalMass) {
			proposedStar->id = starID;
			starID++;
			stars.push_back(proposedStar); // todo: which id?!
			pickedTotalMass += proposedStar->mass;
		}
		else
			break;
	}
	if (stars.size() > 0 && debug)
		std::cout << "Mean mass: " << pickedTotalMass / stars.size() << std::endl;
	if(debug)
		std::cout << "Proposed mass of disk: " << totalMass << " Sampled mass: " << pickedTotalMass << std::endl;
	return stars;
}

std::vector<Star*> InitialConditions::initialMassBulge(double totalMass, int& starID){
	std::vector<Star*> stars;
	Star* proposedStar;
	double pickedTotalMass = 0;
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
					proposedStar = new Star(0, m);
					break;
				}
			}
			else{
				temp = factor2 / (m * ln10) * pow(m, -exponent2);
				if (disAccept(gen) < temp) {
					proposedStar = new Star(0, m);
					break;
				}
			}
		}
		if (pickedTotalMass + proposedStar->mass / 2 < totalMass) {
			proposedStar->id = starID;
			stars.push_back(proposedStar); // todo: which id?!
			pickedTotalMass += proposedStar->mass;
			starID++;
		}
		else
			break;
	}
	if (stars.size()>0 && debug)
		std::cout << "Mean mass: " << pickedTotalMass / stars.size() << std::endl;
	if (debug)
		std::cout << "Proposed mass of Bulge: " << totalMass << " Sampled mass: " << pickedTotalMass << std::endl;
	return stars;
}

void InitialConditions::plummerSphere(std::vector<Star*>& stars, double totalMass, double scaleParameter, double G){
	std::uniform_real_distribution<> dis(0.0, 0.99);//avoid close to singularity
	for (Star* star : stars) {
		double distance = scaleParameter /sqrt(pow(dis(gen), -2. / 3.) - 1);
		star->position = Vec3D::randomVector(distance);
		plummerVelocity(star, scaleParameter, distance, totalMass, G);
	}
}

void InitialConditions::offsetCluster(std::vector<Star*>& stars, Vec3D& offset){
	for (Star* star : stars) {
		star->position += offset;
	}
}

double InitialConditions::plummerEscapeVelocity(double distance, double structuralLength, double totalMass, double G){
	//return sqrt(2.) * pow(distance * distance + structuralLength, -0.25);
	//https://github.com/bacook17/behalf/blob/master/behalf/initialConditions.py
	return sqrt(2. * G * totalMass / structuralLength) *pow(1.+distance*distance/(structuralLength* structuralLength),-0.25);
}

void InitialConditions::plummerVelocity(Star* star, double structuralLength, double distance, double totalMass, double G){
	std::uniform_real_distribution<> dis(0.0, 1.0);
	//Rejection Technique
	std::uniform_real_distribution<> dis_g(0.0, 0.1);
	double q = dis(gen); // random value in range [0,1]
	double g = dis_g(gen); // random value in range [0,0.1]
	while (g > q* q* pow( (1. - q * q), 3.5)) {
		q = dis(gen);
		g = dis_g(gen);
	}
	double velocity = q * plummerEscapeVelocity(distance, structuralLength, totalMass, G)/1.01;
	star->velocity = Vec3D::randomAngles(velocity);
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