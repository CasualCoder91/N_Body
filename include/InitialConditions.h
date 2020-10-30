/**
 * Handles initialization of mass, position and velocity of Stars.
 *
 * @author Alarich Herzner
 * @version 0.2 15.05.2020
*/

#pragma once

#include<vector>
#include <random>

#include "Star.h"
#include "MWPotential.h"
#include "ProgressBar.h"
#include "WangPotential.h"

extern bool debug;

class InitialConditions{
private:
	/**@brief 1km divided by 1pc*/
	//static double kmInpc;
	/**@brief pseudo-random generator used for sampling various distributions*/
	std::mt19937 gen;
	/**@brief The Potential used during simulation. Initial conditions heavily depend on it.*/
	MWPotential* potential;

public:
	/**@brief Prefered constructor. All member variables are initialized */
	InitialConditions(MWPotential* potential);
	/**
	@brief Creates stars with default member variables (mass, position, velocty, acceleration)
	@param [in,out] firstID ID of the first star in the return vector. Gets incremented with every added star.
	@param nStars Amout of stars generated
	*/
	static std::vector<Star*> initStars(int& firstID, int nStars);
	/**
	@brief Creates field stars in the field of view defined by \p focus, \p viewPoint, \p distance and \p angleOfView.
	@param dx The step size along the line of sight in pc
	@param [in,out] firstID ID of the first star in the return vector. Gets incremented with every added star.
	@param focus one point along the line of sight (0,0,0) would be the center of the MW.
	@param viewPoint the position of the observer.
	@param distance how far the observer can see/the render distance
	@param angleOfView the angle of view in degrees
	*/
	std::vector<Star*> initFieldStars(int& firstID, Vec3D focus, Vec3D viewPoint, double distance, double angleOfView);

	double bulgeStarMass(Vec3D focus, Vec3D viewPoint, double distance, double dx, double angleOfView);
	/**
	@brief Sets mass of stars by inverse transform sampling a Salpeter IMF.
	@param [in,out] stars The mass of these stars will be modified.
	@param minMass The smallest possible mass.
	@param maxMass The largest possible mass.
	@param alpha The exponent used in the IMF.
	*/
	double initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha= -2.35);


	double brokenPowerLaw(std::vector<Star*>& stars, std::vector<double> massLimits, std::vector<double> exponents);
	/**
	@brief Creates stars belonging to the disk inside the given cuboid. 
	Sub-cubes are created, the totall mass within which is calculated, the stars are sampled such that the total sampled mass comes close to the value optained by integration.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	@param tlf Top-Left-Front corner of the cube in [kpc]
	@param brf Bottom-Right-Front corner of the cube in [kpc]
	@param depth Depth of the cube in z direction [kpc]
	@param gridResolution Lenght of one sub-cube [kpc]. This dictates the accuracy.
	*/
	/*std::vector<Star*> initDiskStars(int& starID, Vec3D tlf, Vec3D brf, double depth, double gridResolution = 0.001);*/
	/**
	@brief Creates stars belonging to the disk with mass optained through rejection sampling.
	@param totalMass The sum of stellar masses should be equal to the totalMass. In actuality the sum is a bit larger.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	*/
	std::vector<Star*> diskIMF(double totalMass, int& starID);
	/**
	@brief Creates stars belonging to the bulge/spheroid with mass optained through rejection sampling.
	@param totalMass The sum of stellar masses should be equal to the totalMass. In actuality the sum is a bit larger.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	@see bulgeIMF
	*/
	std::vector<Star*> bulgeIMF(double totalMass, int& starID);

	void sampleWang(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement);
	/**
	@brief Sets positions of stars by rejection sampling the density function of the disc.
	@param [in,out] stars The positions of these stars will be modified.
	@param position One corner (typicaly bottom left front) of the volume the positions lie within.
	@param volumeElement leghts of the sides of the cube (typicaly all positive and equal size) [kpc].
	@see sampleBulgePositions
	*/
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement);
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix* transformationMatrix);
	/**
	@brief Adds velocity (circular with dispersion) at given \p position to \p velocity 
	@param [in,out] velocity [km/s]
	@param [in] position [pc]
	@see sampleDiskVelocities
	*/
	void sampleDiskVelocity(Vec3D& velocity, Vec3D& position);
	/**
	@brief Adds velocity (circular with dispersion) at given star positions to corresponding star velocities.
	@param [in,out] stars
	@see sampleDiskVelocity
	*/
	double sampleDiskVelocities(std::vector<Star*> stars);
	/**
	@brief Sets positions of stars by rejection sampling the density function of the bulge.
	@param [in,out] stars The positions of these stars will be modified.
	@param position One corner (typicaly bottom left front) of the volume the positions lie within.
	@param volumeElement leghts of the sides of the cube (typicaly all positive and equal size) [kpc].
	@see sampleDiskPositions
	*/
	void sampleBulgePositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement);
	void sampleBulgePositions(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double coneR, double distance, Matrix* transformationMatrix);
	/**
	@brief Adds velocity (circular with dispersion) at given \p position to \p velocity
	@param [in,out] velocity [km/s]
	@param [in] position [pc]
	@see sampleBulgeVelocities
	*/
	void sampleBulgeVelocity(Vec3D& velocity, Vec3D& position); //todo:test!
	/**
	@brief Adds velocity (circular with dispersion) at given star positions to corresponding star velocities.
	@param [in,out] stars
	@see sampleBulgeVelocity
	*/
	void sampleBulgeVelocities(std::vector<Star*> stars);
	/**
	@brief Sets position and velocity of \p stars overwriting(!) current values. Used to setup clusters.
	@param [in,out] stars Obtain this parameter with for instance initStars(int& firstID)
	@param totalMass The sum of star masses.
	@param scaleParameter size of the cluster core [pc]
	@param G gravitational constant
	*/
	void plummerSphere(std::vector<Star*>& stars, double totalMass, double scaleParameter, double G);
	/**
	@brief Adds \p offset to positions of \p stars 
	@param [in,out] stars
	@param offset [pc]
	*/
	void offsetCluster(std::vector<Star*>& stars, static Vec3D& offset) const;

private:
	friend class Test;
	/** @brief Calculates the local escape velocity [km*s^-1] in the Plummer model */
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass, double G);
	/**
	@brief Sets the (isotropic) velocity of the \p star in a plummer sphere.
	@param [in,out] star
	@param scaleParameter size of the cluster core [pc]
	@param distance of the \p star from the center of the sphere [pc]
	@param totalMass [SolarMass] of the cluster
	@param G gravitational constant
	*/
	void plummerVelocity(Star* star, double scaleParameter, double distance, double totalMass, double G);
	/** @brief calculates the value which is closest to zero in the intervall [a,b] */
	static double closestToZero(double a, double b);
	/** @brief calculates the value which is farthest away from zero in the intervall [a,b] */
	static double farthermostFromZero(double a, double b);
	static void setBoundaries(double& min, double& max);
	static void setBoundaries(Vec3D& min, Vec3D& max);
};

//double diskStarMass(Vec3D focus, Vec3D viewPoint, double distance, double dx, double angleOfView);

//double sampleDiskPositionsNew(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double minDist, double maxDist, double maxR, Matrix* transformationMatrix);

//std::vector<Star*> initFieldStars(int& firstID, Vec3D focus, Vec3D viewPoint, double distance, double dx, double angleOfView);

//void sampleBulgePositionsNew(std::vector<Star*> stars, Vec3D coneBoundaryMin, Vec3D coneBoundaryMax, double minDist, double maxDist, double maxR, Matrix* transformationMatrix);
