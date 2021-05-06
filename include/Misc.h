#pragma once

#include "Constants.h"

double aBrightness(double luminosity, double distance);

//implementation of Mass-Luminosity relation
double luminosity(double mass); //mass in [mSun]

// luminosity == absolute magnitude
//returns apparent magnitude (relative to the luminosity of the sun!)
double apparentMagnitude(double absoluteMagnitude, double distance);

//http://hosting.astro.cornell.edu/academics/courses/astro201/mag_absolute.htm
double absoluteMagnitude(double luminosity);