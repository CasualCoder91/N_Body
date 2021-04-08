#include "Projection.h"

const double  Projection::EtoPJ2000_0[3] = { -0.0548761806632, 0.4941094158461, -0.867666116641649 };
const double  Projection::EtoPJ2000_1[3] = { -0.8734369590164, -0.4448300538949, -0.198076 };
const double  Projection::EtoPJ2000_2[3] = { -0.4838351820817, 0.746982, 0.455984 };
const double* const Projection::EtoP[3] = { EtoPJ2000_0,EtoPJ2000_1,EtoPJ2000_2 };

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

void Projection::GCAtoLSR(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	static double z = 0., s = 0., c = 1.;
	positionOut.x = Constants::positionSun.x - positionIn.x; //pc
	positionOut.y = -positionIn.y; //pc
	positionOut.z = positionIn.z - Constants::positionSun.z; //pc
	velocityOut.x = -velocityIn.x; //km/s
	velocityOut.y = Constants::circularVelocitySun - velocityIn.y; //km/s
	velocityOut.z = velocityIn.z; //km/s
	if (Constants::positionSun.z != 0) { // if not 0 -  need to rotate a bit for GC to have positionOut.z=0
		double t;
		if (z != Constants::positionSun.z) {
			z = Constants::positionSun.z; //pc
			t = hypot(Constants::positionSun.z, Constants::positionSun.x); //pc
			s = Constants::positionSun.z / t; //no unit
			c = Constants::positionSun.x / t; //no unit
		}
		t = positionOut.x; //pc
		positionOut.x = c * t - s * positionOut.z; //pc
		positionOut.z = s * t + c * positionOut.z; //pc
		t = velocityOut.x; //km/s
		velocityOut.x = c * t - s * velocityOut.z; //km/s
		velocityOut.z = s * t + c * velocityOut.z; //km/s
	}
}

void Projection::LSRtoGCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut) {
	static double z = 0., s = 0., c = 1.;
	if (Constants::positionSun.z != 0) { // need to rotate a bit for GC to have rv[3][2]=0
		/*vec6 in = rv[3];*/
		Vec3D positionInLocal = positionIn;
		Vec3D velocityInLocal = velocityIn;
		double t;
		if (z != Constants::positionSun.z) {
			z = Constants::positionSun.z;
			t = hypot(Constants::positionSun.z, Constants::positionSun.x);
			s = Constants::positionSun.z / t;
			c = Constants::positionSun.x / t;
		}
		t = positionInLocal.x;
		positionInLocal.x = c * t + s * positionInLocal.z;
		positionInLocal.z = -s * t + c * positionInLocal.z;
		t = velocityInLocal.x;
		velocityInLocal.x = c * t + s * velocityInLocal.z;
		velocityInLocal.z = -s * t + c * velocityInLocal.z;
		positionOut.x = Constants::positionSun.x - positionInLocal.x;
		positionOut.y = -positionInLocal.y;
		positionOut.z = positionInLocal.z + Constants::positionSun.z;
		velocityOut.x = -velocityInLocal.x;
		velocityOut.y = Constants::circularVelocitySun - velocityInLocal.y;
		velocityOut.z = velocityInLocal.z;
	}
	else {
		positionOut.x = Constants::positionSun.x - positionIn.x;
		positionOut.y = -positionIn.y;
		positionOut.z = positionIn.z;
		velocityOut.x = -velocityIn.x;
		velocityOut.y = Constants::circularVelocitySun - velocityIn.y;
		velocityOut.z = velocityIn.z;
	}
}



void Projection::GCAtoGCP(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	double x2y2 = positionIn.x * positionIn.x + positionIn.y * positionIn.y;
	double x2y2z2 = x2y2 + positionIn.z * positionIn.z;
	positionOut.x = sqrt(x2y2z2);
	positionOut.y = atan2(positionIn.y, positionIn.x);
	positionOut.z = asin(positionIn.z / positionOut.x);
	velocityOut.x = (positionIn.x * velocityIn.x + positionIn.y * velocityIn.y + positionIn.z * velocityIn.z) / positionOut.x;
	velocityOut.y = (velocityIn.x * positionIn.y - positionIn.x * velocityIn.y) / x2y2;
	velocityOut.z = (positionIn.z * (positionIn.x * velocityIn.x + positionIn.y * velocityIn.y) - x2y2 * velocityIn.z) / (x2y2z2 * sqrt(x2y2));
}

void Projection::GCPtoGCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	positionOut.x = positionIn.x * cos(positionIn.z) * cos(positionIn.y);
	positionOut.y = positionIn.x * cos(positionIn.z) * sin(positionIn.y);
	positionOut.z = positionIn.x * sin(positionIn.z);
	velocityOut.x = cos(positionIn.z) * cos(positionIn.y) * velocityIn.x + positionIn.x * cos(positionIn.z) * sin(positionIn.y) * velocityIn.y + positionIn.x * sin(positionIn.z) * cos(positionIn.y) * velocityIn.z;
	velocityOut.y = cos(positionIn.z) * sin(positionIn.y) * velocityIn.x - positionIn.x * cos(positionIn.z) * cos(positionIn.y) * velocityIn.y + positionIn.x * sin(positionIn.z) * sin(positionIn.y) * velocityIn.z;
	velocityOut.z = sin(positionIn.z) * velocityIn.x - positionIn.x * cos(positionIn.z) * velocityIn.z;
}

void Projection::LSRtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	positionOut.x = positionIn.x;
	positionOut.y = positionIn.y;
	positionOut.z = positionIn.z;
	velocityOut.x = velocityIn.x - Constants::velocitySun.x;
	velocityOut.y = velocityIn.y - Constants::velocitySun.y;
	velocityOut.z = velocityIn.z - Constants::velocitySun.z;
}

void Projection::HCAtoLSR(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	positionOut.x = positionIn.x;
	positionOut.y = positionIn.y;
	positionOut.z = positionIn.z;
	velocityOut.x = velocityIn.x + Constants::velocitySun.x;
	velocityOut.y = velocityIn.y + Constants::velocitySun.y;
	velocityOut.z = velocityIn.z + Constants::velocitySun.z;
}

void Projection::HCAtoHEQ(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut) {
	//HCA internal: (kpc, kpc/Myr)
	//HEQ internal: (kpc,radian,kpc/Myr,radian/Myr)
	Vec3D positionInLocal = 0.001 * positionIn; //pc -> kpc
	Vec3D velocityInLocal = 0.001 * velocityIn; //km/s -> kpc/Myr || kpc/Myr = 0.001 * pc/yr = 0.001 * (3,086e+13)km/(3,154e+7)s = 978 km/s

	Vec3D h1,h2;
	h1.x = positionInLocal.x * EtoP[0][0] + positionInLocal.y * EtoP[0][1] + positionInLocal.z * EtoP[0][2];
	h1.y = positionInLocal.x * EtoP[1][0] + positionInLocal.y * EtoP[1][1] + positionInLocal.z * EtoP[1][2];
	h1.z = positionInLocal.x * EtoP[2][0] + positionInLocal.y * EtoP[2][1] + positionInLocal.z * EtoP[2][2];
	h2.x = velocityInLocal.x * EtoP[0][0] + velocityInLocal.y * EtoP[0][1] + velocityInLocal.z * EtoP[0][2];
	h2.y = velocityInLocal.x * EtoP[1][0] + velocityInLocal.y * EtoP[1][1] + velocityInLocal.z * EtoP[1][2];
	h2.z = velocityInLocal.x * EtoP[2][0] + velocityInLocal.y * EtoP[2][1] + velocityInLocal.z * EtoP[2][2];

	double R = hypot(h1.x, h1.y);
	positionOut.x = hypot(R, h1.z);
	double
		ca = (R == 0.) ? 1. : h1.x / R,
		sa = (R == 0.) ? 0. : h1.y / R,
		cd = (positionOut.x == 0.) ? 1. : R / positionOut.x,
		sd = (positionOut.x == 0.) ? 0. : h1.z / positionOut.x,
		temp = ca * h2.x + sa * h2.y;
	positionOut.y = (sa < 0.) ? Constants::pi2 - acos(ca) : acos(ca);
	positionOut.z = asin(sd);
	velocityOut.x = cd * temp + sd * h2.z;
	velocityOut.y = (positionOut.x == 0.) ? 0. : (ca * h2.y - sa * h2.x) / positionOut.x;
	velocityOut.z = (positionOut.x == 0.) ? 0. : (cd * h2.z - sd * temp) / positionOut.x;

	positionOut.x = positionOut.x * 1000; //kpc -> pc
	positionOut.y = positionOut.y * Constants::radInArcsec; //rad -> arcsec
	positionOut.z = positionOut.z * Constants::radInArcsec; //rad -> arcsec
	velocityOut.x = velocityOut.x * 1000; //kpc/Myr -> km/s
	velocityOut.y = velocityOut.y * Constants::radmyrInArcsecyr;
	velocityOut.z = velocityOut.z * Constants::radmyrInArcsecyr;
}

void Projection::HEQtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	Vec3D h1, h2;
	double
		ca = cos(positionIn.y),
		sa = sin(positionIn.y),
		cd = cos(positionIn.z),
		sd = sin(positionIn.z),
		R = cd * positionIn.x,
		va = positionIn.x * velocityIn.y,
		vd = positionIn.x * velocityIn.z,
		temp = cd * velocityIn.x - sd * vd;
	h1.x = ca * R;
	h1.y = sa * R;
	h1.z = sd * positionIn.x;
	h2.x = ca * temp - sa * va;
	h2.y = sa * temp + ca * va;
	h2.z = sd * velocityIn.x + cd * vd;
	positionOut.x = h1.x * EtoP[0][0] + h1.y * EtoP[1][0] + h1.z * EtoP[2][0];
	positionOut.y = h1.x * EtoP[0][1] + h1.y * EtoP[1][1] + h1.z * EtoP[2][1];
	positionOut.z = h1.x * EtoP[0][2] + h1.y * EtoP[1][2] + h1.z * EtoP[2][2];
	velocityOut.x = h2.x * EtoP[0][0] + h2.y * EtoP[1][0] + h2.z * EtoP[2][0];
	velocityOut.y = h2.x * EtoP[0][1] + h2.y * EtoP[1][1] + h2.z * EtoP[2][1];
	velocityOut.z = h2.x * EtoP[0][2] + h2.y * EtoP[1][2] + h2.z * EtoP[2][2];
}

void Projection::HCAtoHGP(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	double R = hypot(positionIn.x, positionIn.y);
	positionOut.x = hypot(R, positionIn.z);
	double
		cl = (R == 0.) ? 1. : positionIn.x / R,
		sl = (R == 0.) ? 0. : positionIn.y / R,
		cb = (positionOut.x == 0.) ? 1. : R / positionOut.x,
		sb = (positionOut.x == 0.) ? 0. : positionIn.z / positionOut.x,
		temp = cl * velocityIn.x + sl * velocityIn.y;
	positionOut.y = (sl < 0.) ? Constants::pi2 - acos(cl) : acos(cl);
	positionOut.z = asin(sb);
	velocityOut.x = cb * temp + sb * velocityIn.z;
	velocityOut.y = (positionOut.x == 0.) ? 0. : (cl * velocityIn.y - sl * velocityIn.x) / positionOut.x;
	velocityOut.z = (positionOut.x == 0.) ? 0. : (cb * velocityIn.z - sb * temp) / positionOut.x;
}

void Projection::HGPtoHCA(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
	double
		cl = cos(positionIn.y),
		sl = sin(positionIn.y),
		cb = cos(positionIn.z),
		sb = sin(positionIn.z),
		R = cb * positionIn.x,
		vl = positionIn.x * velocityIn.y,
		vb = positionIn.x * velocityIn.z,
		temp = cb * velocityIn.x - sb * vb;
	positionOut.x = cl * R;
	positionOut.y = sl * R;
	positionOut.z = sb * positionIn.x;
	velocityOut.x = cl * temp - sl * vl;
	velocityOut.y = sl * temp + cl * vl;
	velocityOut.z = sb * velocityIn.x + cb * vb;
}

void Projection::HCAtoHTP(const Vec3D& positionIn, Vec3D& positionOut, const Vec3D& focus)
{
	//Vec3D focusHTP = focus;
	//focusHTP.x = focus.length();
	//focusHTP.y = atan2(focus.y, focus.x);
	//focusHTP.z = asin(focus.z / focusHTP.x);

	Matrix rotationM = Matrix::rotation(focus, Vec3D(1, 0, 0));
	Vec3D pHCArotated = rotationM * positionIn; //rotate Position such that x points towards focus

	positionOut.x = positionIn.length(); //radius in pc

	positionOut.y = atan2(pHCArotated.y, pHCArotated.x); //declination in rad (atan2 = angle between vector and z axis)
	if (positionOut.y > Constants::pi)
		positionOut.y -= Constants::pi2;
	else if (positionOut.y < -Constants::pi)
	{
		positionOut.y += Constants::pi2;
	}

	positionOut.z = asin(pHCArotated.z / positionOut.x); //ascension in rad
	if(positionOut.z > Constants::pi)
		positionOut.z -= Constants::pi2;
	else if (positionOut.z < -Constants::pi)
	{
		positionOut.z += Constants::pi2;
	}

	positionOut.y = positionOut.y * Constants::radInArcsec; //rad -> arcsec
	positionOut.z = positionOut.z * Constants::radInArcsec; //rad -> arcsec
}

//void Projection::HGPtoHEQ(const Vec3D& positionIn, const Vec3D& velocityIn, Vec3D& positionOut, Vec3D& velocityOut){
//	double dNGP = 0.47347728280415174; //27°7'41.z''  | equatorial coordinates of the north Galactic pole
//	double aNGP = 3.366033268750004;   //12h51m26.28s | equatorial coordinates of the north Galactic pole
//	double lNCP = 2.1630214485816124;  //123°55'55.2''| north celestial pole in Galactic coordinates
//
//	double sd = sin(dNGP)*sin(positionIn.z) + cos(dNGP)*cos(positionIn.z)*cos(lNCP - positionIn.y); //sinus of declination
//	double d = asin(sd);
//	double a = asin(cos(positionIn.z) * sin(lNCP - positionIn.y) / cos(d)) + aNGP;
//
//	positionOut.x = positionIn.x;
//	positionOut.y = a * Constants::radInArcsec;
//	positionOut.z = d * Constants::radInArcsec;
//	return;
//}
