#pragma once
#include "Vec3D.h"
#include <iostream>

class Matrix {
private:
	double m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44;

public:
	Matrix();
	Matrix(double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44);

	static Matrix transformation(Vec3D rotation, Vec3D translation);

	Vec3D operator * (const Vec3D& rhs);

	friend std::ostream& operator<<(std::ostream& output, const Matrix& m) {
		output << m.m11 << " " << m.m12 << " " << m.m13 << std::endl;
		output << m.m21 << " " << m.m22 << " " << m.m23 << std::endl;
		output << m.m31 << " " << m.m32 << " " << m.m33 << std::endl;

		return output;
	}

};