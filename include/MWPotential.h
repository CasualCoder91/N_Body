/**
 * The Milky Way Potential. A combination of the potentials for the black hole, bulge, disk and dark matter halo.
 *
 * @author Alarich Herzner
 * @version 0.3 15.05.2020
 * @todo: remove inheritance from Parameters.
*/

#pragma once

#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include "Vec3D.h"
#include "Star.h"
#include "InOut.h"
#include "LookupTable.h"
#include "ProgressBar.h"
#include "Matrix.h"

#include "Potential/Potential.h"
#include "Potential/Hernquist.h"

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

extern bool debug;

class MWPotential
{
private:
	static const double mMassBlackHole; // SolarMassUnit
	static const double mMassDisk; // SolarMassUnit
	static const double mMassBulge; // SolarMassUnit
	static const double mMassSmallBulge; // SolarMassUnit
	static const double rHalo; // kpc
	static const double mMassHalo; // SolarMassUnit
	static const double mDensityHalo; // SolarMassUnit * pc^-3
	/** @brief 1 km divided by 1 pc */
	//static const double kmInpc;

	static const double velocityDispersionScaleLength;

	static const std::string lookupTableLocation;

public:
	static const double characteristicVelocityBulge; // km/s
	static const double aBulge; // SolarMassUnit
	//static const double aSmallBulge; // SolarMassUnit
	static const double aDisk; // pc
	static const double bDisk; // pc
	static const std::string velocityDistributionBulgeTableFilename;
	LookupTable velocityDistributionBulgeTable;
	static Hernquist bulgePotential;

private:
	static double closestToZero(double a, double b);
	static double gslDensity(double x[], size_t dim, void* p);
	static double gslDensity(double z, void* p);
	static double gslVelocityBulge(double r, void* p);
	static double gslDensityDisk(double x[], size_t dim, void* p);
	static double gslDensityDisk(double z, void* p);

	//Density in cylinder volume along z Axis
	static double gslDensityDiskx(double x, void* p);
	static double gslDensityDisky(double y, void* p);
	static double gslDensityDiskz(double z, void* p);

	/**
	@brief Calculates (numerical integration) the sperical averaged value of the disc potential derived by r.
	This function is used to generate the bulge velocity dispersion
	@param r The radius [pc]
	*/
	static double sphericalAveragedDisc(double r);
	double angularVelocity(double R); // return in s^-1, new Bulge/Halo
	double surfaceDensity(double R); // return in SolarMassUnit*pc^-2, new Bulge
public:
	MWPotential();
	/** @brief the total (all potentials considered) circular velocity [km/s] at the \p position [pc] */
	double circularVelocity(const Vec3D* position);
	/** @brief the circular velocity [km/s] at the \p position [pc] only considereing the disc potential*/
	double circularVelocityDisk(Vec3D* position);
	/** @brief the circular velocity [km/s] at the \p position [pc] only considereing the black hole potential*/
	double circularVelocityBlackHole(Vec3D* position);
	/** @brief the circular velocity [km/s] at the \p position [pc] only considereing the dark matter halo potential*/
	double circularVelocityHalo(Vec3D* position);
	/** @brief the local escape velocity at the given \p position*/
	double escapeVelocity(Vec3D* position);
	/** @brief the local density of the disc*/
	static double densityDisk(double R, double z);
	/** @brief the local density of the disc*/
	static double densityDisk(double x, double y, double z);

	static double denistyWang(double x, double y, double z);
	/**
	@brief Caluclate the surface mass density of the disk at a given radial distance \p R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	double surfaceDensityDisk(double R);
	/**@brief: Mass of the disc inside the \p volumeElement relative to the given \p position*/
	double massDisk(Vec3D position, Vec3D volumeElement);

	double massDisk(Matrix* transformation, double distance, double r);

	//all Potentials
	/**@brief the radial (R) velocity dispersion [km/s] for stars belonging to the disc*/
	double radialVelocityDispersionDisk(double R, double z);
	/**@brief the vertical (z) velocity dispersion [km/s] for stars belonging to the disc*/
	double verticalVelocityDispersion(double R);
	/**@brief the azimuthal velocity dispersion [km/s] for stars belonging to the disc (tangential to circle around the center)*/
	double azimuthalVelocityDispersion(double R, double z);
	/**@brief the average azimuthal velocity [km/s] for stars belonging to the disc (tangential to circle around the center)*/
	double azimuthalStreamingVelocity(Vec3D position);
	/** 
	@brief: the velocity dispersion due to the bulge potential at radius \p r
	@note: This function is resource intensive and should only be used to create a lookup table.
	@see: velocityDistributionBulgeTableValue
	*/
	double velocityDispersionBulge(double r);
	/** @brief: the velocity dispersion due to the bulge potential at radius \p r read from a lookup table */
	double velocityDistributionBulgeTableValue(double r);
	/** @brief: creates a lookup table for the velocity distribution of the bulge. It does **not** initialize the lookuptable for the potential. */
	void generateVelocityDistributionBulgeLookupTable(double rMax);
	/** @brief: adds acceleration based on the potential to the acceleration of the stars. */
	void applyForce(Star* star);
	void applyForce(const Vec3D position, Vec3D& acceleration);

	static double potentialEnergy(const Vec3D& position);
	static double potentialEnergy(double R, double z);


	//returns [km/pc * 1/s]
	double epicyclicFrequency(double R, double z); // in s^-1 new Bulge/Halo

	//static double radialVelocityDispersionBulge(double R, double z); // return in km/s
	//static double infiniteDistributionFunctionBulge(double q); //q = sqrt(-E*characteristicVelocity^-2)
	//static double distributionFunctionBulge(double e); // e = relative Energy = -E
	//static double particleEnergy(Star* star); // Energy per unit mass (not relative energy)
	//static double particleEnergy(Vec3D& position, Vec3D& velocity);
	//static double particleEnergy(Vec3D& position, double velocity);

	//moved into InitialConditions
	//static Vec3D sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
	//static Vec3D sampleBuldge(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

};

/*
	static double gslDensityDiskxNew(double x, void* p);
	static double gslDensityDiskyNew(double y, void* p);
	static double gslDensityDiskzNew(double z, void* p);
	double massDisk(Matrix* transformation, double minDist, double maxDist, double maxR);
*/