#pragma once

#include "SimulationData.h"
#include "Vec3D.h"

# define M_PI 3.14159265358979323846

class Potential : SimulationData
{
private:
	const double massBlackHole = 4e6; // SolarMassUnit
	const double massDisk = 10e10; // SolarMassUnit
	const double aDist = 6.5; //kpc
	const double bDist = 0.26; //kpc
	const double massBuldge = 3.45e10; // SolarMassUnit
	const double aBuldge = 0.7; // SolarMassUnit
	const double rHalo = 16; // kpc
	const double densityHalo = 7e6; // SolarMassUnit*kpc^-3
	const double G = 4.302e-6; // 
	Vec3D position;

public:
	Potential(Vec3D position);
	double circularVelocity(Vec3D* position);
	double circularVelocityDisk(Vec3D* position);
	double circularVelocityBlackHole(Vec3D* position);
	double circularVelocityBuldge(Vec3D* position);
	double circularVelocityHalo(Vec3D* position);


};