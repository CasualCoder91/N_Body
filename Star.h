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

	std::string Dump();//returns all stored data as string
	void Reset(); // set all variables to 0
};

