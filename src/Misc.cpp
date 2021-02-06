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
