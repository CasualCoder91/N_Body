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
	static const double aBulge; // SolarMassUnit
	static const double rHalo; // kpc
	static const double densityHalo; // SolarMassUnit*kpc^-3
	static const double G; // parsec*solarMassUnit^-1*km^2/s^2
	/** @brief 1 km divided by 1 pc */
	static const double kmInpc;

private:
	static double closestToZero(double a, double b);
public:

	static double circularVelocity(Vec3D* position); // return in km/s
	static double circularVelocityDisk(Vec3D* position);
	static double circularVelocityBlackHole(Vec3D* position);
	static double circularVelocityBulge(Vec3D* position);
	static double circularVelocityHalo(Vec3D* position);

	static double densityDisk(double R, double z);
	static double densityDisk(double x, double y, double z);
	/**
	@brief Caluclate the surface mass density of the disk at a given radial distance R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	static double surfaceDensityDisk(double R);
	static Vec3D sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	static double densityBulge(double r);
	static double densityBulge(double x, double y, double z);
	/**
	@brief Caluclate the surface mass density of the bulge at a given radial distance R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	static double surfaceDensityBulge(double R);
	static Vec3D sampleBuldge(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

	//Mass of Disk and Bulge inside volume
	static double frequencyDistribution(Vec3D position, Vec3D volumeElement);
	static double massDisk(Vec3D position, Vec3D volumeElement);
	static double massBulge(Vec3D position, Vec3D volumeElement);

	//all Potentials
	static double angularVelocity(double R); // return in s^-1
	static double surfaceDensity(double R); // return in SolarMassUnit*pc^-2, tested
	static double epicyclicFrequency(double R); // in s^-1
	static double radialVelocityDispersion(double R); // return in km/s
	static double verticalVelocityDispersion(double R); // return in km/s
	static double azimuthalVelocityDispersion(double R); // return in km/s
	static double azimuthalStreamingVelocity(Vec3D position);

	static void applyForce(Star* star);

};