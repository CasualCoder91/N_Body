/**
 * Star containing information about mass, position, velocity and acceleration.
 *
 * @author Alarich Herzner
 * @version 0.9 02.03.2020
*/

#pragma once
#include "Vec3D.h"
class Star{
public:
	Star();
	Star(double mass);
	Star(double mass, Vec3D position);
	Star(double mass, double x, double y, double z);
	double mass;
	Vec3D position;
	Vec3D velocity;
	Vec3D acceleration;

	std::string dump();//returns all stored data as string
	void reset(); // set all variables to 0
};

