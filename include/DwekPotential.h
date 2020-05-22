#pragma once

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_legendre.h>

#include "Vec3D.h"

class DwekPotential {
	double density0 = 0.8;
	double q = 0.6;
	double z0 = 0.4e3; //pc
	double x0 = 1.49e3; //pc
	double y0 = 0.58e3; //pc
	double rs = 1e3; //pc


	double denisty(double x, double y, double z);

	double PotentialNLM(unsigned int n, unsigned int l, unsigned int m, Vec3D position); //Basis function for Potential

};