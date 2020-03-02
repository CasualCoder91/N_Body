#pragma once

#include<vector>
#include <random>

#include "Star.h"
#include "Parameters.h"

class InitialConditions: Parameters {
public:
	static double InitialMass(std::vector<Star*> &stars,int n_Stars=0);
	static void PlummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter
private:
	static double PlummerEscapeVelocity(double distance, double structuralLength, double totalMass);
	static void PlummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
};

