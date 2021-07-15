#pragma once

#include "Constants.h"
#include <gsl/gsl_interp.h>

double aBrightness(double luminosity, double distance);

//implementation of Mass-Luminosity relation
double luminosity(double mass); //mass in [mSun]


//returns apparent magnitude 
double apparentMagnitude(double mass, double distance);

//http://hosting.astro.cornell.edu/academics/courses/astro201/mag_absolute.htm
double absoluteMagnitudeV(double luminosity);

double absoluteMagnitude(double mass);