#include "Projection.h"

//double Projection::kmInpc = 3.086e-13;

Vec3D Projection::projectPosition(const Vec3D& target, const Vec3D& lookAt, const Vec3D& origin, const double fovAngle) {
	Vec3D returnValue = target - origin;
	double length = returnValue.length();
	returnValue = returnValue / length;
	double angleXY = atan2(lookAt.y, lookAt.x) - atan2(returnValue.y, returnValue.x);
	double angleZ_XY = asin(lookAt.z) - asin(returnValue.z);
	returnValue.x = length * tan(angleXY);
	returnValue.y = length * tan(angleZ_XY);
	return returnValue;
}

//http://www.cns.nyu.edu/~david/handouts/motion.pdf
//https://stackoverflow.com/questions/13832505/world-space-to-camera-space
void Projection::project(Vec3D& position, Vec3D& velocity, Vec3D& lookAt, Vec3D& origin) {

	static double focalLength = 100;
	lookAt = lookAt.normalize();

	Vec3D uP = Vec3D(lookAt.y, -lookAt.x, 0);
	Vec3D r = Vec3D::crossProduct(&lookAt, &uP).normalize();
	Vec3D u = Vec3D::crossProduct(&r, &lookAt).normalize();

	Matrix m = Matrix(r.x, r.y, r.z, -r.x * origin.x - r.y * origin.y - r.z * origin.z,
		u.x, u.y, u.z, -u.x * origin.x - u.y * origin.y - u.z * origin.z,
		-lookAt.x, -lookAt.y, -lookAt.z, lookAt.x * origin.x + lookAt.y * origin.y + lookAt.z * origin.z,
		0, 0, 0, 1);

	Vec3D positionCam = m * position;

	position = Vec3D(focalLength * positionCam.x / positionCam.z, focalLength * positionCam.y / positionCam.z, 0);
	velocity = focalLength / positionCam.z * velocity + focalLength * velocity.z / (positionCam.z * positionCam.z) * positionCam;

	//Matrix mT = Matrix(1, 0, 0, -origin.x, 0, 1, 0, -origin.y, 0, 0, 1, -origin.z, 0, 0, 0, 1);
	//Matrix m = mR * mT;

	return;
}

Vec3D Projection::projectVelocity(const Vec3D& target, const Vec3D& lookAt, const double fovAngle) {
	Vec3D returnValue = target;
	double length = returnValue.length();
	returnValue = returnValue / length;
	double angleXY = atan2(lookAt.y, lookAt.x) - atan2(returnValue.y, returnValue.x);
	double angleZ_XY = asin(lookAt.z) - asin(returnValue.z);
	returnValue.x = length * tan(angleXY);
	returnValue.y = length * tan(angleZ_XY);
	return returnValue;
}