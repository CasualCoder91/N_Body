/**
 * The Milky Way Potential. A combination of the potentials for the black hole, bulge, disk and dark matter halo.
 *
 * @author Alarich Herzner
 * @version 0.2 16.04.2020
*/

#pragma once

#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include "SimulationData.h"
#include "Vec3D.h"
#include "Star.h"

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

class Potential : SimulationData
{
private:
	static const double mMassBlackHole; // SolarMassUnit
	static const double mMassDisk; // SolarMassUnit
	static const double aDisk; //kpc
	static const double bDisk; //kpc
	static const double mMassBulge; // SolarMassUnit
	//static const double mMassBulge2; // SolarMassUnit
	static const double aBulge; // SolarMassUnit
	//static const double aBulge2; // SolarMassUnit
	static const double rHalo; // kpc
	//static const double densityHalo; // SolarMassUnit*kpc^-3
	static const double mMassHalo; // SolarMassUnit
	static const double G; // parsec*solarMassUnit^-1*km^2/s^2
	/** @brief 1 km divided by 1 pc */
	static const double kmInpc;

	static const double velocityDispersionScaleLength;
public:
	static const double characteristicVelocityBulge; // km/s

private:
	static double closestToZero(double a, double b);
	//averaged Disk potential at radius r
	
public:
	static double sphericalAveragedDisk(double r);
	static double circularVelocity(Vec3D* position); // return in km/s
	static double circularVelocityDisk(Vec3D* position);
	static double circularVelocityBlackHole(Vec3D* position);
	static double circularVelocityBulge(Vec3D* position); // with new Bulge!
	static double circularVelocityHalo(Vec3D* position);

	static double escapeVelocity(Vec3D* position);

	static double densityDisk(double R, double z);
	static double densityDisk(double x, double y, double z);
	/**
	@brief Caluclate the surface mass density of the disk at a given radial distance R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	static double surfaceDensityDisk(double R);
	//static Vec3D sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	static double densityBulge(double r); //new Bulge
	static double densityBulge(double R, double z); //new Bulge
	static double densityBulge(double x, double y, double z); //new Bulge
	/**
	@brief Caluclate the surface mass density of the bulge at a given radial distance R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	static double surfaceDensityBulge(double R); //new Bulge
	//static Vec3D sampleBuldge(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	//Mass of Disk and Bulge inside volume
	static double frequencyDistribution(Vec3D position, Vec3D volumeElement);
	static double massDisk(Vec3D position, Vec3D volumeElement);
	static double massBulge(Vec3D position, Vec3D volumeElement);

	//all Potentials
	static double angularVelocity(double R); // return in s^-1, new Bulge/Halo
	static double surfaceDensity(double R); // return in SolarMassUnit*pc^-2, new Bulge
	static double epicyclicFrequency(double R, double z); // in s^-1 new Bulge/Halo
	static double radialVelocityDispersionDisk(double R, double z); // return in km/s
	static double verticalVelocityDispersion(double R); // return in km/s
	static double azimuthalVelocityDispersion(double R, double z); // return in km/s
	static double azimuthalStreamingVelocity(Vec3D position);

	//todo: test this, check Dimensions!
	static double velocityDistributionBulge(double r); // working with new potential?
	static double meanVelocityBulge(double r);
	//static double radialVelocityDispersionBulge(double R, double z); // return in km/s
	//static double infiniteDistributionFunctionBulge(double q); //q = sqrt(-E*characteristicVelocity^-2)
	//static double distributionFunctionBulge(double e); // e = relative Energy = -E
	//static double particleEnergy(Star* star); // Energy per unit mass (not relative energy)
	//static double particleEnergy(Vec3D& position, Vec3D& velocity);
	//static double particleEnergy(Vec3D& position, double velocity);
	static void applyForce(Star* star);
	//static double potentialEnergy(Vec3D& position);
	//static double potentialEnergy(double R, double z);

};