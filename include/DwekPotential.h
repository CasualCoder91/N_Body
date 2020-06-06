#pragma once

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_monte_plain.h>
#include <math.h>       /* pow */

#include "Vec3D.h"

class DwekPotential {
	static const double density0;
	static const double q;
	static const double z0; //pc
	static const double x0; //pc
	static const double y0; //pc
	static const double rs; //pc
	static const double mass;

	static double density(double x, double y, double z);
	static double gslDensity(double x[], size_t dim, void* p);
	//x[0] = (r-1/r+1), x[1] = cos(theta), x[2] = phi
	static double expansionCoeffIntegrand(double x[], size_t dim, void* p);

	double PotentialNLM(unsigned int n, unsigned int l, unsigned int m, Vec3D position); //Basis function for Potential
public:
	static double distributionFunction(Vec3D position, Vec3D velocity); //do: angle!
	static double ANLM(unsigned int n, unsigned int l, unsigned int m);
	static double totalMass(double min, double max);
};