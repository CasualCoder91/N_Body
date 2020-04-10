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

class InitialConditions {
private:
	double G;
	/** @brief Amount of stars in the cluster */
	int nStars;

public:
	InitialConditions(SimulationData* parameters);
	/**
	@brief Creates cluster stars with default member variables (mass, position, velocty, acceleration)
	@param firstID The ID given to the first star in the vector. Subsequent stars get higher IDs
	*/
	std::vector<Star*> initStars(int firstID);
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
	@param tlf Top-Left-Front corner of the cube in [kpc]
	@param tlf Bottom-Right-Front corner of the cube in [kpc]
	@param depth Depth of the cube in z direction [kpc]
	@param gridResolution Lenght of one sub-cube [kpc]. This dictates the accuracy.
	*/
	std::vector<Star*> initDiskStars(int firstID, Vec3D tlf, Vec3D brf, double depth, Potential* potential, double gridResolution = 0.001);
	/**
	@brief Creates stars belonging to the disk with mass optained through rejection sampling.
	@param totalMass The sum of stellar masses should be equal to the totalMass. In actuality the sum is a bit larger.
	*/
	std::vector<Star*> massDisk(double totalMass);
	/**
	@brief Sets positions of stars by rejection sampling the density function of the disc.
	@param [in,out] stars The possitions of these stars will be modified.
	@param position One corner (typicaly bottom left front) of the volume the positions lie within.
	@param volumeElement leghts of the sides of the cube (typicaly all positive and equal size) [kpc].
	@param potential Pointer of the Potential defining the density distribution of the disk.
	@see sampleBulgePositions
	*/
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement, Potential* potential);
	/**
	@brief Sets positions of stars by rejection sampling the density function of the bulge.
	@param [in,out] stars The possitions of these stars will be modified.
	@param position One corner (typicaly bottom left front) of the volume the positions lie within.
	@param volumeElement leghts of the sides of the cube (typicaly all positive and equal size) [kpc].
	@param potential Pointer of the Potential defining the density distribution of the disk.
	@see sampleDiskPositions
	*/
	double sampleBulgePositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement, Potential* potential);

	void plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter

	void setNStars(int N);
private:
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass);
	void plummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
	double closestToZero(double a, double b);
};

