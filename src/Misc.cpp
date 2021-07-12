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

double apparentMagnitude(double luminosity, double distance) {
	//https://en.wikipedia.org/wiki/Absolute_magnitude
	return absoluteMagnitude(luminosity) - 5.0 + 5.0 * log10(distance);
}

double absoluteMagnitude(double luminosity){
	//parameter luminosity is relativ to luminosity of the sun (=3.828e26W)
	//division by fixed luminosity (=3.0128e28W) leads to factor in log10
	//https://en.wikipedia.org/wiki/Luminosity#Relationship_to_magnitude
	return - 2.5*log10(luminosity* 0.01270578863515666);
}
