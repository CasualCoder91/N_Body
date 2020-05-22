#include "..\include\DwekPotential.h"

double DwekPotential::denisty(double x, double y, double z){
	double r1 = pow(gsl_pow_2(gsl_pow_2(x / x0) + gsl_pow_2(y / y0)) + gsl_pow_4(z / z0), 0.25);
	double r2 = sqrt((gsl_pow_2(q) * (gsl_pow_2(x) + gsl_pow_2(y)) + gsl_pow_2(z)) / gsl_pow_2(z0));
	return density0 * (exp(-gsl_pow_2(r1) / 2) + pow(r2, -1.85) * exp(-r2));
}

double DwekPotential::PotentialNLM(unsigned int n, unsigned int l, unsigned int m, Vec3D position){
	double s = position.length() / rs;
	return pow(s, l) / pow(1 + s, 2 * l + 1) * gsl_sf_gegenpoly_n(n, 2 * (double)l + 1.5, (s - 1) / (s + 1)) * gsl_sf_legendre_sphPlm(l, m, cos(position.theta())) * cos(m * position.phi());
}
