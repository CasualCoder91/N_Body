#include "Matrix.h"

Matrix::Matrix(){
	this->m11 = 0;
	this->m12 = 0;
	this->m13 = 0;
	this->m14 = 0;
	this->m21 = 0;
	this->m22 = 0;
	this->m23 = 0;
	this->m24 = 0;
	this->m31 = 0;
	this->m32 = 0;
	this->m33 = 0;
	this->m34 = 0;
	this->m41 = 0;
	this->m42 = 0;
	this->m43 = 0;
	this->m44 = 0;
}

Matrix::Matrix(double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33, double m34, double m41, double m42, double m43, double m44){
	this->m11 = m11; this->m12 = m12; this->m13 = m13; this->m14 = m14;
	this->m21 = m21; this->m22 = m22; this->m23 = m23; this->m24 = m24;
	this->m31 = m31; this->m32 = m32; this->m33 = m33; this->m34 = m34;
	this->m41 = m41; this->m42 = m42; this->m43 = m43; this->m44 = m44;
}

Matrix Matrix::transformation(Vec3D rotation, Vec3D translation){
	rotation = rotation.normalize();
	Vec3D zAxis = Vec3D(0, 0, 1);
	double angle = acos(zAxis * rotation)*0.5;
	Vec3D b = Vec3D::crossProduct(&zAxis, &rotation);
	b = b.normalize();
	double q0 = cos(angle), q1 = sin(angle) * b.x, q2 = sin(angle)*b.y, q3 = sin(angle)*b.z;

	double m11 = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
	double m12 = 2 * (q1 * q2 - q0 * q3);
	double m13 = 2 * (q1 * q3 + q0 * q2);
	double m14 = translation.x;
	double m21 = 2 * (q2 * q1 + q0 * q3);
	double m22 = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
	double m23 = 2 * (q2 * q3 - q0 * q1);
	double m24 = translation.y;
	double m31 = 2 * (q3 * q1 - q0 * q2);
	double m32 = 2 * (q3 * q2 + q0 * q1);
	double m33 = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
	double m34 = translation.z;
	double m41 = 0; 
	double m42 = 0; 
	double m43 = 0;
	double m44 = 1;
	return Matrix(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44);
}

Vec3D Matrix::operator*(const Vec3D& rhs){
	return Vec3D(m11*rhs.x+m12*rhs.y+m13*rhs.z+m14, 
		m21 * rhs.x + m22 * rhs.y + m23 * rhs.z + m24,
		m31 * rhs.x + m32 * rhs.y + m33 * rhs.z + m34);
}
