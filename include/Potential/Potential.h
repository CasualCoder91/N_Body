#pragma once

class Potential {

public:
	virtual double density(double r) = 0;
	virtual double density(double R, double z) = 0;
	virtual double density(double x, double y, double z) = 0;
};