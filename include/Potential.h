#pragma once

#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_math.h>
#include "SimulationData.h"
#include "Vec3D.h"

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

class Potential : SimulationData
{
private:
	const double massBlackHole = 4e6; // SolarMassUnit
	const double mMassDisk = 10e10; // SolarMassUnit
	const double aDisk = 6.5; //kpc
	const double bDisk = 0.26; //kpc
	const double mMassBulge = 3.45e10; // SolarMassUnit
	const double aBulge = 0.7; // SolarMassUnit
	const double rHalo = 16; // kpc
	const double densityHalo = 7e6; // SolarMassUnit*kpc^-3
	const double G = 4.302e-6; // 
	Vec3D position;
	//Vec3D volumeElement = Vec3D(0.1, 0.1, 0.1); // dx,dy,dz in kpc

private:
	double closestToZero(double a, double b);
public:
	Potential(Vec3D position);
	double circularVelocity(Vec3D* position);
	double circularVelocityDisk(Vec3D* position);
	double circularVelocityBlackHole(Vec3D* position);
	double circularVelocityBulge(Vec3D* position);
	double circularVelocityHalo(Vec3D* position);

	double densityDisk(double R, double z);
	double densityDisk(double x, double y, double z);
	Vec3D sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	double densityBulge(double r);
	double densityBulge(double x, double y, double z);
	Vec3D sampleBuldge(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	//Mass of Disk and Bulge inside volume
	double frequencyDistribution(Vec3D position, Vec3D volumeElement);
	double massDisk(Vec3D position, Vec3D volumeElement);
	double massBulge(Vec3D position, Vec3D volumeElement);
};