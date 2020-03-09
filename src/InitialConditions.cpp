#include "InitialConditions.h"

double InitialConditions::initialMass(std::vector<Star*> &stars,int n_Stars)
{
	double mass = 0;
	double totalMass = 0;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 0.1);
	if (n_Stars > 0) {
		for (int i = 0; i < n_Stars; i++) {
			mass = dis(gen);
			mass = 0.08 + (0.19 * pow(mass, 1.55) + 0.05 * pow(mass, 0.6)) / pow(1 - mass, 0.58);
			Star* star = new Star(stars.size(),mass);
			stars.push_back(star);
			totalMass += mass;
		}
	}
	else {
		for (Star* star : stars) {
			mass = dis(gen);
			mass = 0.08 + (0.19 * pow(mass, 1.55) + 0.05 * pow(mass, 0.6)) / pow(1 - mass, 0.58);
			star->mass = mass;
			totalMass += mass;
		}
	}
	return totalMass;
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

double InitialConditions::plummerEscapeVelocity(double distance, double structuralLength, double totalMass){
	//return sqrt(2.) * pow(distance * distance + structuralLength, -0.25);
	//https://github.com/bacook17/behalf/blob/master/behalf/initialConditions.py
	return sqrt(2. * Parameters::G * totalMass / structuralLength) *pow(1.+distance*distance/(structuralLength* structuralLength),-0.25);
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
