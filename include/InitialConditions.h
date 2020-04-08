#pragma once

#include<vector>
#include <random>

#include "Star.h"
#include "Parameters.h"
#include "Potential.h"

class InitialConditions {
private:
	double G;
	int nStars;

public:
	InitialConditions(SimulationData* parameters);
	std::vector<Star*> initStars(int firstID);
	double initialMass(std::vector<Star*> &stars);
	double initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha= -2.35);
	std::vector<Star*> initDiskStars(int firstID, Vec3D tlf, Vec3D brf, double depth, Potential* potential);
	std::vector<Star*> massDisk(double totalMass); //rejection sampling
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement, Potential* potential);
	void plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter
private:
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass);
	void plummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
	double rangeZero(double a, double b);
};

