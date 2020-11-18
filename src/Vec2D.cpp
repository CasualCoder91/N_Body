#include "Vec2D.h"

Vec2D::Vec2D() {
	x = y = 0;
}

Vec2D::Vec2D(double x, double y) {
	this->x = x;
	this->y = y;
}

Vec2D::~Vec2D() {}

double Vec2D::length() const{
	return sqrt(this->x * this->x + this->y * this->y);
}
