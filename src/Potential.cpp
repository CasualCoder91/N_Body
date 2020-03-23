#include "..\include\Potential.h"

Potential::Potential(Vec3D position){
	this->position = position;
}



double Potential::circularVelocity(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y,2);
	double velocity = sqrt(G * massBlackHole * R2 / pow(R2 + pow(z, 2), 1.5)
		+ G * massDisk * R2 / pow(pow(aDist + sqrt(pow(bDist, 2) + pow(z, 2)), 2) + R2, 1.5)
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
	double velocity = sqrt( G * massDisk * R2 / pow(pow(aDist + sqrt(pow(bDist, 2) + pow(z, 2)), 2) + R2, 1.5));
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
