/**
 * Handles initialization of mass, position and velocity of Stars.
 *
 * @author Alarich Herzner
 * @version 0.1 09.04.2020
*/

#pragma once

#include<vector>
#include <random>

#include "Star.h"
#include "Parameters.h"
#include "Potential.h"
#include "ProgressBar.h"

extern bool debug;

class InitialConditions : Parameters {
private:
	/**@brief 1km divided by 1pc*/
	static double kmInpc;
	/**@brief pseudo-random generator used for sampling various distributions*/
	std::mt19937 gen;
	/**@brief The Potential used during simulation. Initial conditions heavily depend on it.*/
	Potential* potential;

public:
	/**@brief Prefered constructor. All member variables are initialized */
	InitialConditions(Parameters* parameters, Potential* potential);
	/**
	@brief Creates stars with default member variables (mass, position, velocty, acceleration)
	@param [in,out] firstID ID of the first star in the return vector. Gets incremented with every added star.
	*/
	std::vector<Star*> initStars(int& firstID);
	/**
	@brief Creates stars with default member variables (mass, position, velocty, acceleration)
	@param [in,out] firstID ID of the first star in the return vector. Gets incremented with every added star.
	@param nStars Amout of stars generated
	@see initStars(int& firstID)
	*/
	static std::vector<Star*> initStars(int firstID, int nStars);
	/**
	@brief Creates stars with default member variables (mass, position, velocty, acceleration)
	@param [in,out] firstID ID of the first star in the return vector. Gets incremented with every added star.
	@param Amount of stars given by Parameters::nStars
	@see initStars(int& firstID)
	*/
	std::vector<Star*> initFieldStars(int& firstID);
	/**
	@brief Sets mass of stars by inverse transform sampling a Salpeter IMF.
	@param [in,out] stars The mass of these stars will be modified.
	@param minMass The smallest possible mass.
	@param maxMass The largest possible mass.
	@param alpha The exponent used in the IMF.
	*/
	double initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha= -2.35);
	/**
	@brief Creates stars belonging to the disk inside the given cuboid. 
	Sub-cubes are created, the totall mass within which is calculated, the stars are sampled such that the total sampled mass comes close to the value optained by integration.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	@param tlf Top-Left-Front corner of the cube in [kpc]
	@param brf Bottom-Right-Front corner of the cube in [kpc]
	@param depth Depth of the cube in z direction [kpc]
	@param gridResolution Lenght of one sub-cube [kpc]. This dictates the accuracy.
	*/
	std::vector<Star*> initDiskStars(int& starID, Vec3D tlf, Vec3D brf, double depth, double gridResolution = 0.001);
	/**
	@brief Creates stars belonging to the disk with mass optained through rejection sampling.
	@param totalMass The sum of stellar masses should be equal to the totalMass. In actuality the sum is a bit larger.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	*/
	std::vector<Star*> massDisk(double totalMass, int& starID);
	/**
	@brief Creates stars belonging to the bulge/spheroid with mass optained through rejection sampling.
	@param totalMass The sum of stellar masses should be equal to the totalMass. In actuality the sum is a bit larger.
	@param [in,out] starID ID of the first star in the return vector. Gets incremented with every added star.
	@see initialMassBulge
	*/
	std::vector<Star*> initialMassBulge(double totalMass, int& starID);
	/**
	@brief Sets positions of stars by rejection sampling the density function of the disc.
	@param [in,out] stars The positions of these stars will be modified.
	@param position One corner (typicaly bottom left front) of the volume the positions lie within.
	@param volumeElement leghts of the sides of the cube (typicaly all positive and equal size) [kpc].
	@see sampleBulgePositions
	*/
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement);
	/**
	@brief Adds velocity (circular with dispersion) at given /rev position to /ref velocity 
	@param [in,out] velocity [km/s]
	@param [in] position [pc]
	@see sampleDiskVelocities
	*/
	void sampleDiskVelocity(Vec3D& velocity, Vec3D& position);
	/**
	@brief Adds velocity (circular with dispersion) at given star positions to corresponding star velocities.
	@param [in,out] velocity [km/s]
	@param [in] position [pc]
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

	void sampleBulgeVelocity(Vec3D& velocity, Vec3D& position); //todo:test!
	void sampleBulgeVelocities(std::vector<Star*> stars);

	void plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter

	void offsetCluster(std::vector<Star*>& stars, Vec3D& offset);

	void setNStars(int N);
private:
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass); // in km*s^-1
	void plummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
	static double closestToZero(double a, double b);
	static double farthermostFromZero(double a, double b);
};

