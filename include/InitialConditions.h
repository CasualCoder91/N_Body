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
	@paran minMass The smallest possible mass.
	@paran maxMass The largest possible mass.
	@paran alpha The exponent used in the IMF.
	*/
	double initialMassSalpeter(std::vector<Star*>& stars, double minMass, double maxMass, double alpha= -2.35);
	/**
	@brief Creates stars belonging to the disk inside the given cuboid. 
	Sub-cubes are created, the totall mass within which is calculated, the stars are sampled such that the total sampled mass comes close to the value optained by integration.
	@paran tlf Top-Left-Front corner of the cube in [kpc]
	@paran tlf Bottom-Right-Front corner of the cube in [kpc]
	@paran depth Depth of the cube in z direction [kpc]
	*/
	std::vector<Star*> initDiskStars(int firstID, Vec3D tlf, Vec3D brf, double depth, Potential* potential);
	/**
	@brief Creates stars belonging to the disk with mass optained through rejection sampling.
	@param totalMass Mass which the summ of star masses should be equal to.
	*/
	std::vector<Star*> massDisk(double totalMass); //rejection sampling
	double sampleDiskPositions(std::vector<Star*> stars, Vec3D position, Vec3D volumeElement, Potential* potential);
	void plummerSphere(std::vector<Star*>& stars, double structuralLength, double totalMass); // structuralLength = a = softening parameter

	void setNStars(int N);
private:
	double plummerEscapeVelocity(double distance, double structuralLength, double totalMass);
	void plummerVelocity(Star* star, double structuralLength, double distance, double totalMass);
	double rangeZero(double a, double b);
};

