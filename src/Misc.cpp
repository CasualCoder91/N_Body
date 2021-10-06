#include "..\include\Misc.h"

double aBrightness(double luminosity, double distance){
	return luminosity/(Constants::pi4*pow(distance,2));
}

double luminosity(double mass){
	if (mass < 0.43) {
		return 0.23 * pow(mass, 2.3);
	}
	else if (mass < 2.) {
		return pow(mass, 4);
	}
	else if (mass < 55) {
		return 1.4 * pow(mass, 3.5);
	}
	else{
		return 32000. * mass;
	}
}

double apparentMagnitude(double mass, double distance) {
	//https://en.wikipedia.org/wiki/Absolute_magnitude
	return absoluteMagnitude(mass) - 5.0 + 5.0 * log10(distance);
}

double absoluteMagnitudeV(double luminosity){
	//parameter luminosity is relativ to luminosity of the sun 
	//https://en.wikipedia.org/wiki/Luminosity#Relationship_to_magnitude
	return 4.74 - 2.5*log10(luminosity);
}

double absoluteMagnitude(double mass){
	//Eric Mamajek Version 2021.03.02
	//https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
	if (mass > 19.8) {
		double lum = luminosity(mass);
		return absoluteMagnitudeV(lum);
	}
	static const size_t n = 81;
	static const double Msun[] = { 0.001,0.075,0.076,0.077,0.078,0.079,0.08,0.085,0.088,0.09,0.093,0.102,0.123,0.162,0.184,0.23,0.27,0.37,0.4,0.44,0.47,0.5,0.54,0.57,0.59,0.62,0.64,0.69,0.7,0.73,0.78,0.82,0.86,0.88,0.9,0.94,0.95,0.97,0.98,0.985,0.99,1,1.03,1.06,1.08,1.13,1.18,1.21,1.25,1.33,1.38,1.44,1.46,1.5,1.61,1.75,1.77,1.81,1.83,1.86,1.88,1.93,1.98,2.05,2.18,2.68,2.75,3.38,3.92,4.3,4.7,5.1,5.4,6,7.3,10,11,15,17.7,18.5,19.8 };
	static const double M_Ks[] = { 20,11,10.77,10.55,10.5,10.4,10.3,9.92,9.81,9.7,9.5,9.22,8.8,8.2,7.93,7.36,7.1,6.55,6.18,5.98,5.75,5.64,5.36,5.15,5.01,4.95,4.81,4.56,4.397,4.247,4.1,3.938,3.889,3.827,3.693,3.532,3.495,3.409,3.345,3.319,3.282,3.236,3.12,3.043,2.947,2.915,2.76,2.579,2.5,2.291,2.188,2.116,2.045,1.941,1.836,1.756,1.702,1.694,1.655,1.447,1.607,1.587,1.172,1.07,0.949,0.648,0.621,0.254,-0.075,-0.192,-0.433,-0.553,-0.708,-0.956,-1.198,-1.848,-2.126,-2.587,-2.942,-3.073,-3.2 };
	gsl_interp* workspace = gsl_interp_alloc(gsl_interp_linear, n);
	gsl_interp_init(workspace, Msun, M_Ks, n);
	double mag = gsl_interp_eval(workspace, Msun, M_Ks, mass, NULL);	 
	gsl_interp_free(workspace);
	return mag;
}
