#include "..\include\Potential.h"


double rangeZero(double a, double b) {
	double temp = 0;
	if ((a < 0 && b> 0) || (a > 0 && b<0))
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

Potential::Potential(Vec3D position){
	this->position = position;
}



double Potential::circularVelocity(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y,2);
	double velocity = sqrt(G * massBlackHole * R2 / pow(R2 + pow(z, 2), 1.5)
		+ G * massDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5)
		+ G * massBuldge * R2 / (r * pow(aBuldge + r, 2))
		+ 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 3) * log(r / rHalo + 1) / pow(R2+pow(z,2), 1.5)
		- 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 2) / (pow(r, 2) * (r / rHalo + 1)));
	return velocity;
}

double Potential::circularVelocityDisk(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt( G * massDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5));
	return velocity;
}

double Potential::circularVelocityBlackHole(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * massBlackHole * R2 / pow(R2 + pow(z, 2), 1.5));
	return velocity;
}

double Potential::circularVelocityBuldge(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * massBuldge * R2 / (r * pow(aBuldge + r, 2)));
	return velocity;
}

double Potential::circularVelocityHalo(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(4 * M_PI * G * densityHalo * R2 * pow(rHalo, 3) * log(sqrt(R2 + pow(z, 2)) / rHalo + 1) / pow(R2 + pow(z, 2), 1.5)
		- 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 2) / ((R2+pow(z, 2)) * (sqrt(R2 + pow(z, 2)) / rHalo + 1)));
	return velocity;
}

double Potential::densityDisk(double R, double z){
	double z2b2 = pow(z, 2) + pow(bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	return pow(bDisk, 2)* massDisk / (4 * M_PI) * (aDisk * pow(R, 2) + (aDisk + 3 * sz2b2 * pow(aDisk + sz2b2, 2))) / (pow(pow(R, 2) + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
}

double Potential::densityDisk(double x, double y, double z){
	double R = sqrt(pow(x, 2) + pow(y, 2));
	return densityDisk(R,z);
}

Vec3D Potential::sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
	double smallestx = rangeZero(xMin, xMax);
	double smallesty = rangeZero(yMin, yMax);
	double smallestz = rangeZero(zMin, zMax);

	double k = densityDisk(smallestx, smallesty, smallestz);
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(xMin, xMax);
	std::uniform_real_distribution<> disy(yMin, yMax);
	std::uniform_real_distribution<> disz(zMin, zMax);
	std::uniform_real_distribution<> disaccept(0, k);

	while (true) {
		double x = disx(gen);
		double y = disy(gen);
		double z = disz(gen);

		double accept = disaccept(gen);
		double temp = densityDisk(x, y, z);

		if (accept < temp)
			return Vec3D(x, y, z);
	}
	return Vec3D();
}
