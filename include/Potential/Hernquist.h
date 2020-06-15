#pragma once

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_plain.h>
#include <stdio.h> //fprintf
#include <iostream>

#include "Vec3D.h"
#include "Star.h"

class Hernquist {

	static double mMass;
	static double mScaleLength; // SolarMassUnit
	static double G;

public:
	Hernquist(double mass, double scaleLength, double G);

	/** @brief the local density*/
	static double density(double r);
	static double density(double R, double z);
	static double density(double x, double y, double z);
	static double gslDensity(double x[], size_t dim, void* p);
	static double gslDensity(double z, void* p);

	double potential(double r);

	/** @brief part of gradient used in MWPotential to apply force to a star*/
	double forceTemp(double r);

	/**@brief: Mass inside the \p volumeElement relative to the given \p position*/
	double mass(Vec3D position, Vec3D volumeElement);
	/**
	@brief Caluclate the surface mass density at a given radial distance R.
	The GSL function gsl_integration_qagiu is used.
	@param R The radius [pc] at which the surface mass density is calculated.
	@return The calculated surface mass density [SolarMassUnit*pc^-2].
	*/
	double surfaceDensity(double R);

	/** @brief the circular velocity [km/s] at the \p position [pc]*/
	double circularVelocity(Vec3D* position);

	double escapeVelocity(double r);
	double escapeVelocity(Vec3D* position);

	/**
	@brief Potential derived from r with G=1. Used for velocity dispersion.
	*/
	double potentialdr(double r);
	/**
	@brief Potential 2x derived from r with G=1. Used for epicyclic frequency.
	*/
	double potentialdr2(double r);
};