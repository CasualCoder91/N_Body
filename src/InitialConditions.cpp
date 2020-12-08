#include "InitialConditions.h"



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

void InitialConditions::setBoundaries(double& min, double& max){
	if ((min < 0 && max> 0) || (min > 0 && max < 0)) {
		if (abs(min) > abs(max)) {
			max = min;
		}
		min = 0;
	}
	else if (abs(min) > abs(max)) {
		double temp = min;
		min = max;
		max = temp;
	}
}

void InitialConditions::setBoundaries(Vec3D& min, Vec3D& max){
	setBoundaries(min.x, max.x);
	setBoundaries(min.y, max.y);
	setBoundaries(min.z, max.z);
	//if (min.x == 0 && min.y == 0 && min.z == 0) //remove singularity in bulge density instead
	//	min.x = 1;
}

InitialConditions::InitialConditions(MWPotential* potential):gen((std::random_device())()){
	this->potential = potential;
}

std::vector<Star*> InitialConditions::initStars(int& firstID,int nStars){
	std::vector<Star*> stars = {};
	int max = firstID + nStars;
	for (; firstID < max; firstID++) {
		Star* star =  new Star(firstID);
		stars.push_back(star);
	}
	return stars;
}

//Cone: https://www.math24.net/calculation-volumes-triple-integrals/
std::vector<Star*> InitialConditions::initFieldStars(int& starID, const Vec3D& focus, const Vec3D& viewPoint, double distance, double angleOfView){
	std::cout << "Initializing field stars" << std::endl;
	Vec3D direction = (focus - viewPoint).normalize();

	angleOfView = angleOfView * Constants::degInRad; //convert degrees in rad
	double coneR = distance * tan(0.5*angleOfView);
	Matrix transformationMatrix = Matrix::transformation(direction, viewPoint);
	Vec3D coneBoundaryMin = transformationMatrix * Vec3D(-coneR, -coneR, 0);
	Vec3D coneBoundaryMax = transformationMatrix * Vec3D(coneR, coneR, distance * 1.10); //1.01 to make sure boundary is not inside cone.
	setBoundaries(coneBoundaryMin, coneBoundaryMax);

	std::vector<Star*> fieldStars; //return vector

	double diskMass = potential->massDisk(&transformationMatrix, distance, coneR);
	std::vector<Star*> diskStars = diskIMF(diskMass, starID);
	if (diskStars.size() > 0) {
		sampleDiskPositions(diskStars, coneBoundaryMin, coneBoundaryMax, coneR, distance, &transformationMatrix); //test
		sampleDiskVelocities(diskStars);
		fieldStars.insert(std::end(fieldStars), std::begin(diskStars), std::end(diskStars));
	}
	double bulgeMass = potential->bulgePotential.mass(&transformationMatrix, distance, coneR);
	std::vector<Star*> bulgeStars = bulgeIMF(bulgeMass, starID);
	if (bulgeStars.size() > 0) {
		sampleBulgePositions(bulgeStars, coneBoundaryMin, coneBoundaryMax, coneR, distance, &transformationMatrix); //test
		sampleBulgeVelocities(bulgeStars);
		fieldStars.insert(std::end(fieldStars), std::begin(bulgeStars), std::end(bulgeStars));
	}
	return fieldStars;
}


double InitialConditions::bulgeStarMass(const Vec3D& focus, const Vec3D& viewPoint, double distance, double dx, double angleOfView){
	std::cout << "Calculating Bulge Star Mass" << std::endl;
	Vec3D direction = (focus - viewPoint).normalize();
	int nSteps = static_cast<int>((distance + dx) / dx);
	ProgressBar progressBar = ProgressBar(0, static_cast<float>(nSteps), true);
	double totalMass = 0;
	//#pragma omp parallel for
	for (int step = 1; step <= nSteps; ++step) {//steps along direction (line of sight)
		double r = 0.5 * step * dx * tan(angleOfView); //distance from line of sight at current step
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
		totalMass += potential->bulgePotential.mass(corner, volumeElement);
		if (volumeElement.length() < potential->aBulge) {
			totalMass += potential->massDisk(corner, volumeElement);
		}
		progressBar.Update(static_cast<float>(step));
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

double InitialConditions::brokenPowerLaw(std::vector<Star*>& stars, const std::vector<double>& massLimits, const std::vector<double>& exponents){
	size_t nIntervals = massLimits.size() - 1;
	double totalMass = 0;
	for (double exponent : exponents) {
		if (exponent == 1) {
			std::cout << "Exponentvalue 1 in brokenPowerLaw not possible" << std::endl;
			return -1;
		}
	}
	std::vector<double> continuity(nIntervals); //k[i]
	continuity[0] = pow(massLimits[1], exponents[0]);
	continuity[1] = pow(massLimits[1], exponents[1]);
	for (size_t i = 2; i < nIntervals; ++i) {
		continuity[i] = continuity[i - 1] * pow(massLimits[i], exponents[i] - exponents[i - 1]);
	}

	std::vector<double> exponentTemp(nIntervals); // 1 - alpha 
	double constant = 0; // normalization constant A
	for (size_t i = 0; i < nIntervals; i++) {
		exponentTemp[i] = 1-exponents[i];
		constant += (continuity[i]*(pow(massLimits[i + 1], exponentTemp[i]) - pow(massLimits[i], exponentTemp[i]))) / exponentTemp[i];
	}
	constant = 1 / constant;

	std::vector<double>inverseTemps = { 0 };
	double inverseTemp = 0;
	for (size_t i = 0; i < nIntervals; i++) {
		inverseTemp += constant* continuity[i] / exponentTemp[i] * (pow(massLimits[i + 1], exponentTemp[i]) - pow(massLimits[i], exponentTemp[i]));
		inverseTemps.push_back(inverseTemp);
	}

	std::uniform_real_distribution<> dis(0.0, 1.0);
	for (Star* star : stars) {
		double sample = dis(gen); // sample = y
		for (size_t i = 0; i < nIntervals; i++) {
			if (sample < inverseTemps[i + 1]) {
				star->mass = pow((sample - inverseTemps[i]) * exponentTemp[i] / (constant * continuity[i]) + pow(massLimits[i], exponentTemp[i]), 1 / exponentTemp[i]);
				totalMass += star->mass;
				break;
			}
		}
	}
	return totalMass;
}

void InitialConditions::sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement) {
	//get Limits for "Accept" distribution. Depends on Volume. Largest density is found at positions closest to 0
	double smallestx = closestToZero(position.x, position.x + volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);
	double largestx = farthermostFromZero(position.x, position.x + volumeElement.x);
	double largesty = farthermostFromZero(position.y, position.y + volumeElement.y);
	double largestz = farthermostFromZero(position.z, position.z + volumeElement.z);
	double acceptUpperLimit = potential->densityDisk(smallestx, smallesty, smallestz);
	double acceptLowerLimit = potential->densityDisk(largestx, largesty, largestz);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);

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
	return;
}

void InitialConditions::sampleDiskPositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix* transformationMatrix) {
	double acceptUpperLimit = potential->densityDisk(coneBoundaryMin.x, coneBoundaryMin.y, coneBoundaryMin.z);
	double acceptLowerLimit = potential->densityDisk(coneBoundaryMax.x, coneBoundaryMax.y, coneBoundaryMax.z);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);

	std::uniform_real_distribution<> disx(-coneR, coneR);
	std::uniform_real_distribution<> disy(-coneR, coneR);
	std::uniform_real_distribution<> disz(0, distance);

	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);

			//double yBound = sqrt(coneR * coneR - x * x);
			//std::uniform_real_distribution<> disy(-yBound, yBound);
			//double y = disy(gen);
			//std::uniform_real_distribution<> disz(distance / coneR * sqrt(x * x + y * y), distance);
			//double z = disz(gen);

			double y = disy(gen);
			double z = disz(gen);
			double RTest = sqrt(x * x + y * y);

			if (RTest< coneR && z>distance / coneR * RTest) {
				Vec3D trialPosition = Vec3D(x, y, z);
				trialPosition = *transformationMatrix * trialPosition;

				double accept = disaccept(gen);
				double temp = potential->densityDisk(trialPosition.x, trialPosition.y, trialPosition.z);
				if (temp > acceptUpperLimit || temp < acceptLowerLimit) {
					std::cerr << "Uhoh issue in sampleDiskPositions" << std::endl;
					std::cerr << "coneBoundaryMin: " << coneBoundaryMin.print() << std::endl;
					std::cerr << "coneBoundaryMax: " << coneBoundaryMax.print() << std::endl;
					std::cerr << "trialPosition: " << trialPosition.print() << std::endl;
					//std::cin.clear();
					//std::cin.get();
				}
				if (accept < temp && trialPosition.length()>1) {
					star->position = Vec3D(trialPosition.x, trialPosition.y, trialPosition.z);
					break;
				}
			}
		}
	}
	return;
}

void InitialConditions::sampleDiskVelocity(Vec3D& velocity, Vec3D& position){

	double R = gsl_hypot(position.x,position.y);

	std::normal_distribution<> zVelocityDistribution{ 0,potential->verticalVelocityDispersion(R) };
	double vz = zVelocityDistribution(gen);

	double rDispersion = potential->radialVelocityDispersionDisk(R,position.z);
	double vR = 0;
	if (rDispersion > 0) {
		std::normal_distribution<> radialVelocityDistribution{ 0,rDispersion };
		vR = radialVelocityDistribution(gen);
	}
	else if(debug)
		std::cout << "critical error in InitialConditions::sampleDiskVelocity()\n";

	double aDispersion = potential->azimuthalVelocityDispersion(R, position.z);
	double va = 0;
	if (aDispersion > 0) {
		std::normal_distribution<> azimuthalVelocityDistribution{ potential->azimuthalStreamingVelocity(position),aDispersion };
		va = azimuthalVelocityDistribution(gen);
	}
	else {
		if (debug)
			std::cout << "critical error in InitialConditions::sampleDiskVelocity()\n";
		va = potential->azimuthalStreamingVelocity(position);
	}

	double theta = position.theta();
	Vec3D sampledVelocity = Vec3D(vR * sin(theta) + va * sin(theta + M_PI_2), vR * cos(theta) + va * cos(theta + M_PI_2), vz);
	double escapeVelocity = potential->escapeVelocity(&position);
	if (escapeVelocity < sampledVelocity.length()) {
		if(debug)
			std::cout << "sampleDiskVelocity: sampled velocity larger than escape velocity" << std::endl;
		sampledVelocity = Vec3D::randomAngles(escapeVelocity);
	}
	velocity += sampledVelocity;
}

void InitialConditions::sampleDiskVelocities(std::vector<Star*> stars){
	for (Star* star : stars) {
		sampleDiskVelocity(star->velocity, star->position);
	}
	return;
}

void InitialConditions::sampleBulgePositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement) {
	double smallestx = closestToZero(position.x, position.x + volumeElement.x);
	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
	double smallestz = closestToZero(position.z, position.z + volumeElement.z);
	double largestx = farthermostFromZero(position.x, position.x + volumeElement.x);
	double largesty = farthermostFromZero(position.y, position.y + volumeElement.y);
	double largestz = farthermostFromZero(position.z, position.z + volumeElement.z);
	if (smallestx == 0 && smallesty == 0 && smallestz == 0) {
		smallestx = 0.1;
	}
	double acceptUpperLimit = potential->bulgePotential.density(smallestx, smallesty, smallestz);
	double acceptLowerLimit = potential->bulgePotential.density(largestx, largesty, largestz);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);


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

	//const double lamda = 1 / (MWPotential::aBulge * 2);
	//std::exponential_distribution<double> disr(lamda);
	//std::exponential_distribution<double> disaccept(lamda);
	for (Star* star : stars) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);
			double accept = disaccept(gen);
			double temp = potential->bulgePotential.density(x, y, z);

			if (accept < temp) {
				star->position = Vec3D(x, y, z);
				break;
			}
		}
	}
	return; //todo: return average velocity maybe?
}

void InitialConditions::sampleBulgePositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix* transformationMatrix){

	double acceptUpperLimit = potential->bulgePotential.density(coneBoundaryMin.x, coneBoundaryMin.y, coneBoundaryMin.z);
	double acceptLowerLimit = potential->bulgePotential.density(coneBoundaryMax.x, coneBoundaryMax.y, coneBoundaryMax.z);
	//create distribution with calculated limits
	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);

	std::uniform_real_distribution<> disx(-coneR, coneR);
	std::uniform_real_distribution<> disy(-coneR, coneR);
	std::uniform_real_distribution<> disz(0, distance);

	#pragma omp parallel for
	for(int i = 0; i < stars.size(); ++i) {
		while (true) {
			double x = disx(gen);
			double y = disy(gen);
			double z = disz(gen);
			double RTest = sqrt(x * x + y * y);
			if (RTest< coneR && z>distance / coneR * RTest) {
				Vec3D trialPosition = Vec3D(x, y, z);
				trialPosition = *transformationMatrix * trialPosition;

				double accept = disaccept(gen);
				double temp = potential->bulgePotential.density(trialPosition.x, trialPosition.y, trialPosition.z);
				if (temp > acceptUpperLimit || temp < acceptLowerLimit) {
					std::cout << "Uhoh issue in sampleBulgePositions" << std::endl;
					std::cout << "coneBoundaryMin: " << coneBoundaryMin.print() << std::endl;
					std::cout << "coneBoundaryMax: " << coneBoundaryMax.print() << std::endl;
					std::cout << "trialPosition: " << trialPosition.print() << std::endl;
					//	std::cin.clear();
					//	std::cin.get();
				}
				if (accept < temp) {
					stars[i]->position = Vec3D(trialPosition.x, trialPosition.y, trialPosition.z);
					break;
				}
			}
		}
	}
	return; //todo: return average velocity maybe?
}

void InitialConditions::sampleBulgeVelocity(Vec3D& velocity, Vec3D& position){
	double delta = potential->velocityDistributionBulgeTableValue(position.length());
	/*double vCirc = potential->circularVelocity(&position);*/
	double patternSpeed = 0.06; //km*s-1*pc-1
	double maxSpeed = 400;
	//std::normal_distribution<> velocityDistribution{ 0,delta };
	//double vRand = velocityDistribution(gen);

	std::uniform_real_distribution<> disaccept(0, 1);
	std::uniform_real_distribution<> vd(-maxSpeed, maxSpeed);
	//std::normal_distribution<double> rDist{ 0,delta };
	//velocity += Vec3D::randomVector(rDist(gen));
	//velocity.x += patternSpeed * position.y;
	//velocity.y -= patternSpeed * position.x;


	//std::uniform_real_distribution<> vyd(-maxSpeed + patternSpeed * position.y, maxSpeed - patternSpeed * position.x);
	//std::uniform_real_distribution<> vzd(-maxSpeed, maxSpeed);

	while (true) {
		double vx = vd(gen);
		double vy = vd(gen);
		double vz = vd(gen);
		double v2 = vx*vx+ vy*vy +vz*vz;
		double d2 = gsl_pow_2(delta);
		double accept = disaccept(gen);
		double temp = 4 * M_PI * pow(1 / (2 * M_PI * d2), 1.5) * v2 * exp(-v2 / (2 * d2));

		if (temp > 1) {
			std::cout << "sample Error" << std::endl;
		}

		if (accept < temp) {
			velocity = Vec3D(vx, vy, vz);
			double theta = position.theta();
			velocity.x -= patternSpeed * position.y * sin(theta);
			velocity.y += patternSpeed * position.x * cos(theta);
			//Vec3D meanVelocity = Vec3D(patternSpeed * position.y, -patternSpeed * position.x,0);
			//velocity = (velocity.length() - meanVelocity.length()) * velocity.normalize();
			//velocity = velocity + meanVelocity;

			double escapeVel = potential->escapeVelocity(&velocity);
			if (velocity.length() < escapeVel)
				break;
			else if (debug) {
				std::cout << "sampleBulgeVelocity(): star too fast" << std::endl;
			}
		}
	}







	//old junk
	//double  phi = position.theta();
	//double  theta = position.phi();

	//velocity = Vec3D(vCirc * cos(theta + M_PI_2) + cos(phi)*cos(theta)*vRand, vCirc * sin(theta + M_PI_2) + cos(phi) * sin(theta) * vRand, sin(phi)*vRand);
	//velocity = Vec3D(cos(phi)*cos(theta)*vRand + patternSpeed * position.y, cos(phi) * sin(theta) * vRand -patternSpeed*position.x, sin(phi)*vRand);
}

void InitialConditions::sampleBulgeVelocities(std::vector<Star*> stars){
	for (Star* star : stars) {
		sampleBulgeVelocity(star->velocity, star->position);
	}
}

std::vector<Star*> InitialConditions::diskIMF(double totalMass, int& starID){
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

std::vector<Star*> InitialConditions::bulgeIMF(double totalMass, int& starID){
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

void InitialConditions::offsetCluster(std::vector<Star*>& stars, const Vec3D& offset) const{
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
	star->velocity = Vec3D::randomVector(velocity);
}


//void InitialConditions::sampleBulgeVelocity(Vec3D& velocity, Vec3D& position){
//	//std::random_device rd;  //Will be used to obtain a seed for the random number engine
//	//std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//	//double escapeVelocity = sqrt(-2 * MWPotential::potentialEnergy(position));
//	//std::uniform_real_distribution<> distributionVelocity(0.0, escapeVelocity); // Use escape velocity instead!
//	//double maximumEnergy = MWPotential::particleEnergy(position, Vec3D(0, 0, 0));
//	//double maxDF = MWPotential::distributionFunctionBulge(-maximumEnergy);//(gsl_pow_2(MWPotential::characteristicVelocityBulge)-1); // -1 to avoid division by zero
//	//std::uniform_real_distribution<> distributionDF(0.0, maxDF);
//	//while (true) {
//	//	double mVelocity = distributionVelocity(gen);
//	//	double randomDF = distributionDF(gen);
//	//	double energy = MWPotential::particleEnergy(position, mVelocity);
//	//	double attemptDF = MWPotential::distributionFunctionBulge(-energy);
//	//	if (randomDF < attemptDF) {
//	//		position += Vec3D::randomVector(mVelocity);
//	//		return;
//	//	}
//	//}
//	double R = gsl_hypot(position.x, position.y);
//	//std::cout << position.print() << std::endl;
//	std::normal_distribution<> velocityDistribution{ 0,MWPotential::radialVelocityDispersionBulge(R, position.z) };
//	double vz = velocityDistribution(gen);
//	double vR = velocityDistribution(gen);
//	//std::normal_distribution<> velocityDistribution2{ 200,2*MWPotential::radialVelocityDispersionBulge(R, position.z) };
//	double va = velocityDistribution(gen);
//
//	double theta = atan2(position.y, position.x);
//
//	velocity += Vec3D(vR * cos(theta) + va * cos(theta + M_PI_2), vR * sin(theta) + va * sin(theta + M_PI_2), vz);
//}

//double InitialConditions::diskStarMass(Vec3D focus, Vec3D viewPoint, double distance, double dx, double angleOfView) {
//	std::cout << "Calculating Disk Star Mass" << std::endl;
//	Vec3D direction = (focus - viewPoint).normalize();
//	int nSteps = (distance + dx) / dx;
//	ProgressBar progressBar = ProgressBar(0, nSteps, true);
//	double totalMass = 0;
//	//#pragma omp parallel for
//	for (int step = 1; step <= nSteps; ++step) {//steps along direction (line of sight)
//		double r = 0.5 * step * dx * tan(angleOfView); //distance from line of sight at current step
//		double aBoid = 2 * r;
//		if (aBoid < 1) //if cubes are smaller 1pc^3 density is aproximated 0
//			continue;
//		Vec3D rVec = sqrt(2) * r * Vec3D::crossProduct(&Vec3D(-1, 1, -1), &direction).normalize();
//		rVec.x = -abs(rVec.x);
//		rVec.y = -abs(rVec.y);
//		rVec.z = -abs(rVec.z);
//		Vec3D corner = viewPoint + direction * ((double)step - 1) * dx + rVec;
//		//std::cout << corner.print() << std::endl;
//		Vec3D volumeElement = direction * dx - 2 * rVec;
//		totalMass += potential->massDisk(corner, volumeElement);
//		progressBar.Update(step);
//		progressBar.Print();
//	}
//	return totalMass;
//}

//double InitialConditions::sampleDiskPositionsNew(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double minDist, double maxDist, double maxR, Matrix* transformationMatrix) {
//	double acceptUpperLimit = potential->densityDisk(coneBoundaryMin.x, coneBoundaryMin.y, coneBoundaryMin.z);
//	double acceptLowerLimit = potential->densityDisk(coneBoundaryMax.x, coneBoundaryMax.y, coneBoundaryMax.z);
//	//create distribution with calculated limits
//	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);//lower limit is new -> speedup
//
//	std::uniform_real_distribution<> disz(minDist, maxDist);
//
//	for (Star* star : stars) {
//		while (true) {
//			double z = disz(gen);
//			double xBound = z * maxR / maxDist;
//			std::uniform_real_distribution<> disx(-xBound, xBound);
//			double x = disx(gen);
//			double yBound = sqrt(xBound * xBound - x * x);
//			std::uniform_real_distribution<> disy(-yBound, yBound);
//			double y = disy(gen);
//			Vec3D trialPosition = Vec3D(x, y, z);
//			trialPosition = *transformationMatrix * trialPosition;
//
//			double accept = disaccept(gen);
//			double temp = potential->densityDisk(trialPosition.x, trialPosition.y, trialPosition.z);
//			if (temp > acceptUpperLimit || temp < acceptLowerLimit) {
//				std::cerr << "Uhoh issue in sampleDiskPositionsNew" << std::endl;
//				std::cerr << "coneBoundaryMin: " << coneBoundaryMin.print() << std::endl;
//				std::cerr << "coneBoundaryMax: " << coneBoundaryMax.print() << std::endl;
//				std::cerr << "trialPosition: " << trialPosition.print() << std::endl;
//				std::cin.clear();
//				std::cin.get();
//			}
//			if (accept < temp) {
//				star->position = trialPosition;
//				break;
//			}
//		}
//	}
//	return 0; //todo: return average velocity maybe?
//}

//std::vector<Star*> InitialConditions::initFieldStars(int& starID, Vec3D focus, Vec3D viewPoint, double distance, double dx, double angleOfView) {
//	std::cout << "Initializing field stars" << std::endl;
//	Vec3D direction = (focus - viewPoint).normalize();
//	int nSteps = 1;
//	double minDist = 0;
//	double maxDist = distance;
//	if (dx != 0) {
//		nSteps = (distance + dx) / dx;
//	}
//	else
//		dx = distance;
//	std::vector<Star*> fieldStars; //return vector
//	for (int step = 1; step <= nSteps; ++step) {
//		minDist = (step - 1) * dx;
//		maxDist = step * dx;
//		double minR = minDist * tan(0.5 * angleOfView * Constants::degInRad);
//		double maxR = maxDist * tan(0.5 * angleOfView * Constants::degInRad);
//		Matrix transformationMatrix = Matrix::transformation(direction, viewPoint);
//		Vec3D coneBoundaryMin = transformationMatrix * Vec3D(-maxR, -maxR, minDist * 0.95);
//		Vec3D coneBoundaryMax = transformationMatrix * Vec3D(maxR, maxR, maxDist * 1.05); //factor to make sure boundary is not inside cone.
//		setBoundaries(coneBoundaryMin, coneBoundaryMax);
//
//		double diskMass = potential->massDisk(&transformationMatrix, minDist, maxDist, maxR);
//		std::vector<Star*> diskStars = diskIMF(diskMass, starID);
//		if (diskStars.size() > 0) {
//			sampleDiskPositionsNew(diskStars, coneBoundaryMin, coneBoundaryMax, minDist, maxDist, maxR, &transformationMatrix);
//			//sampleDiskVelocities(diskStars);
//			fieldStars.insert(std::end(fieldStars), std::begin(diskStars), std::end(diskStars));
//		};
//		double bulgeMass = potential->bulgePotential.mass(&transformationMatrix, minDist, maxDist, maxR);
//		std::vector<Star*> bulgeStars = bulgeIMF(bulgeMass, starID);
//		if (bulgeStars.size() > 0) {
//			sampleBulgePositionsNew(bulgeStars, coneBoundaryMin, coneBoundaryMax, minDist, maxDist, maxR, &transformationMatrix);
//			//sampleBulgeVelocities(bulgeStars);
//			fieldStars.insert(std::end(fieldStars), std::begin(bulgeStars), std::end(bulgeStars));
//		};
//	}
//	return fieldStars;
//}

//void InitialConditions::sampleBulgePositionsNew(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double minDist, double maxDist, double maxR, Matrix* transformationMatrix) {
//
//	double acceptUpperLimit = potential->bulgePotential.density(coneBoundaryMin.x, coneBoundaryMin.y, coneBoundaryMin.z);
//	double acceptLowerLimit = potential->bulgePotential.density(coneBoundaryMax.x, coneBoundaryMax.y, coneBoundaryMax.z);
//	//create distribution with calculated limits
//	std::uniform_real_distribution<> disaccept(acceptLowerLimit, acceptUpperLimit);//lower limit is new -> speedup
//
//	std::uniform_real_distribution<> disz(minDist, maxDist);
//
//	for (Star* star : stars) {
//		while (true) {
//			double z = disz(gen);
//			double xBound = z * maxR / maxDist;
//			std::uniform_real_distribution<> disx(-xBound, xBound);
//			double x = disx(gen);
//			double yBound = sqrt(xBound * xBound - x * x);
//			std::uniform_real_distribution<> disy(-yBound, yBound);
//			double y = disy(gen);
//			Vec3D trialPosition = Vec3D(x, y, z);
//			trialPosition = *transformationMatrix * trialPosition;
//
//			double accept = disaccept(gen);
//			double temp = potential->bulgePotential.density(trialPosition.x, trialPosition.y, trialPosition.z);
//			if (temp > acceptUpperLimit || temp < acceptLowerLimit) {
//				std::cout << "Uhoh issue in sampleBulgePositionsNew" << std::endl;
//				std::cout << "coneBoundaryMin: " << coneBoundaryMin.print() << std::endl;
//				std::cout << "coneBoundaryMax: " << coneBoundaryMax.print() << std::endl;
//				std::cout << "trialPosition: " << trialPosition.print() << std::endl;
//				std::cin.clear();
//				std::cin.get();
//			}
//			if (accept < temp) {
//				star->position = trialPosition;
//				break;
//			}
//		}
//	}
//	return; //todo: return average velocity maybe?
//}

//std::vector<Star*> InitialConditions::initDiskStars(int& starID, Vec3D tlf, Vec3D brf, double depth, double gridResolution){
//	std::vector<Star*> stars;
//	if (depth < 0) {
//		throw "initDiskStars: Depth must be >= 0";
//	}
//	for (double z = tlf.z; z < tlf.z + depth; z += gridResolution) {
//		std::cout << "z = " << z << " max z = " << tlf.z + depth << std::endl;
//		for (double x = tlf.x; x < brf.x; x += gridResolution) {
//			for (double y = brf.y; y < tlf.y; y += gridResolution) {
//				Vec3D position = Vec3D(x, y, z);
//				Vec3D volumeElement = Vec3D(gridResolution, gridResolution, gridResolution);
//				double massInCell = potential->massDisk(position, volumeElement);
//				std::vector<Star*> starsInCell = diskIMF(massInCell, starID); //stars with mass
//				sampleDiskPositions(starsInCell, position, volumeElement);
//				stars.insert(stars.end(), starsInCell.begin(), starsInCell.end());
//			}
//		}
//	}
//
//	return stars;
//}

//void InitialConditions::sampleWang(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement) {
//	double maximumVelocity = 500; // replace with escape velocity
//	double patternSpeed = 0.06; //km*s-1*pc-1
//	double smallestx = closestToZero(position.x, position.x + volumeElement.x);
//	double smallesty = closestToZero(position.y, position.y + volumeElement.y);
//	double smallestz = closestToZero(position.z, position.z + volumeElement.z);
//	double largestx = farthermostFromZero(position.x, position.x + volumeElement.x);
//	double largesty = farthermostFromZero(position.y, position.y + volumeElement.y);
//	double largestz = farthermostFromZero(position.z, position.z + volumeElement.z);
//	double minMass = stars.at(0)->mass;
//	double maxMass = stars.at(0)->mass;
//	for (Star* star : stars) {
//		if (star->mass < minMass)
//			minMass = star->mass;
//		if (star->mass > maxMass)
//			maxMass = star->mass;
//	}
//	double acceptUpperLimit = WangPotential::distributionFunction(minMass, Vec3D(smallestx, smallesty, smallestz), Vec3D(0, 0, 0)) * 4; //todo: do it better
//	//double acceptLowerLimit = WangPotential::distributionFunction(maxMass,Vec3D(largestx, largesty, largestz), Vec3D(maximumVelocity + patternSpeed * abs(largesty), maximumVelocity + patternSpeed * abs(largestx), maximumVelocity));
//
//	//std::cout << "r:" << position.length() << " acceptUpperLimit:" << acceptUpperLimit << std::endl;
//
//	if (isnan(acceptUpperLimit))
//		double fail = 1;
//
//	std::uniform_real_distribution<> disaccept(0, acceptUpperLimit);
//
//
//	double x1_low = std::min(position.x, position.x + volumeElement.x);
//	double x1_high = std::max(position.x, position.x + volumeElement.x);
//	double x2_low = std::min(position.y, position.y + volumeElement.y);
//	double x2_high = std::max(position.y, position.y + volumeElement.y);
//	double x3_low = std::min(position.z, position.z + volumeElement.z);
//	double x3_high = std::max(position.z, position.z + volumeElement.z);
//	std::uniform_real_distribution<> disx(x1_low, x1_high); //position
//	std::uniform_real_distribution<> disy(x2_low, x2_high);
//	std::uniform_real_distribution<> disz(x3_low, x3_high);
//
//
//	//double vcl = 250 / (1 + pow(100 / (sqrt(pow(x1_high,2)+pow(x2_high,2))), 0.2));
//	//double theta = atan2(position.y, position.x);
//	//double r = sqrt(pow(position.x, 2) + pow(position.y, 2));
//	//std::uniform_real_distribution<> disvx(-100,100); //velocity
//	//std::uniform_real_distribution<> disvy(vcl - 50,vcl+ 50); //velocity
//	//std::uniform_real_distribution<> disvz(-100, 100); //velocity
//	std::uniform_real_distribution<> disvx(-maximumVelocity, maximumVelocity); //velocity
//	std::uniform_real_distribution<> disvy(-maximumVelocity, maximumVelocity); //velocity
//	std::uniform_real_distribution<> disvz(-maximumVelocity, maximumVelocity); //velocity
//
//	for (Star* star : stars) {
//		while (true) {
//			double x = disx(gen);
//			double y = disy(gen);
//			double z = disz(gen);
//
//			double vx = disvx(gen);
//			double vy = disvy(gen);
//			double vz = disvz(gen);
//			//Vec3D velCyl = Vec3D(vx, vy, vz);
//			//theta = atan2(y, x);
//			//r = sqrt(pow(x, 2) + pow(y, 2));
//			//Vec3D velCart = Vec3D(velCyl.x * cos(theta) - r * sin(theta) * velCyl.y, velCyl.x * sin(theta) + r * cos(theta) * velCyl.y, velCyl.z);
//
//			double accept = disaccept(gen);
//			double temp = WangPotential::distributionFunction(star->mass, Vec3D(x, y, z), Vec3D(vx, vy, vz));
//
//			if (temp > acceptUpperLimit) {
//				std::cout << "sampleWang Error" << std::endl;
//			}
//
//			if (accept < temp) {
//				vx += patternSpeed * position.y;
//				vy -= patternSpeed * position.x;
//				//double angle = -0.349066; //20 degree
//				//double xtemp = vx;
//				//double ytemp = vy;
//				//vx = xtemp * cos(angle) - ytemp * sin(angle);
//				//vy = xtemp * sin(angle) + ytemp * cos(angle);
//				star->position = Vec3D(x, y, z);
//				star->velocity = Vec3D(vx, vy, vz);
//				break;
//			}
//		}
//	}
//	return;
//}

//double InitialConditions::brokenPowerLawOld(std::vector<Star*>& stars, const std::vector<double>& massLimits, const std::vector<double>& exponents) {
//	size_t nIntervals = massLimits.size() - 1;
//	double totalMass = 0;
//	for (double exponent : exponents) {
//		if (exponent == 1) {
//			std::cout << "Exponentvalue 1 in brokenPowerLaw not possible" << std::endl;
//			return -1;
//		}
//	}
//
//	double constant = 0;
//	for (size_t i = 0; i < nIntervals; i++) {
//		double exponentTemp = -exponents[i] + 1;
//		constant += (pow(massLimits[i + 1], exponentTemp) - pow(massLimits[i], exponentTemp)) / exponentTemp;
//	}
//	constant = 1 / constant;
//
//	std::vector<double>inverseTemps = { 0 };
//	double inverseTemp = 0;
//	for (size_t i = 0; i < nIntervals; i++) {
//		double exponentTemp = -exponents[i] + 1;
//		inverseTemp -= constant / exponentTemp * (pow(massLimits[i + 1], exponentTemp) - pow(massLimits[i], exponentTemp));
//		inverseTemps.push_back(inverseTemp);
//	}
//
//	std::uniform_real_distribution<> dis(0.0, 1.0);//avoid close to singularity	
//	for (Star* star : stars) {
//		double sample = dis(gen);
//		for (int i = 0; i < nIntervals; i++) {
//			if (sample < -inverseTemps[i + 1]) {
//				double exponentTemp = -exponents[i] + 1;
//				star->mass = pow((inverseTemps[i] + sample) * exponentTemp / constant + pow(massLimits[i], exponentTemp), 1 / exponentTemp);
//				totalMass += star->mass;
//				break;
//			}
//		}
//	}
//	return totalMass;
//}
