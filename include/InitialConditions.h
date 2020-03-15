#pragma once

#include<vector>
#include <random>

#include "Star.h"
#include "Parameters.h"

class InitialConditions: Parameters {
public:
	std::vector<Star*> initStars(int firstID, int n_Stars);
	double initialMass(std::vector<Star*> &stars);
	double initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha= -2.35);
	void plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter
private:
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass);
	void plummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
};

