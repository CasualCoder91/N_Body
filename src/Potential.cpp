#include "..\include\Potential.h"

const double Potential::mMassBlackHole = 4e6; // SolarMassUnit
const double Potential::mMassDisk = 10e10; // SolarMassUnit
const double Potential::aDisk = 6.5e3; //pc
const double Potential::bDisk = 0.26e3; //pc
const double Potential::mMassBulge = 3.45e10; // SolarMassUnit
const double Potential::aBulge = 0.7e3; //pc
const double Potential::rHalo = 16e3; //pc
const double Potential::densityHalo = 7e-3; // SolarMassUnit*pc^-3
const double Potential::G = 4.302e-3; // parsec*solarMassUnit^-1*km^2/s^2
const double Potential::kmInpc = 3.086e-13;
const double Potential::velocityDispersionScaleLength = aDisk / 4.0;

struct gslDensityDiskParams { double mMassDisk; double aDisk; double bDisk; };
struct gslDensityDiskQagsParams { double mMassDisk; double aDisk; double bDisk; double R; };
struct gslDensityBulgeParams { double mMassBulge; double aBulge; };
struct gslDensityBulgeQagsParams { double mMassBulge; double aBulge; double R; };
struct gslDensityParams { double mMassBulge; double aBulge; double mMassDisk; double aDisk; double bDisk; };
struct gslDensityQagsParams {double mMassBulge; double aBulge; double mMassDisk; double aDisk; double bDisk; double R; };

//what about black hole?!
double gslDensity(double x[], size_t dim, void* p) {

	struct gslDensityParams* fp = (struct gslDensityParams*)p;

	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double z2b2 = pow(x[2], 2) + pow(fp->bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R = gsl_hypot(x[0], x[1]);
	double r = gsl_hypot3(x[0], x[1], x[2]);
	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * pow(R, 2) + (fp->aDisk + 3.0 * sz2b2) * pow(fp->aDisk + sz2b2, 2)) / (pow(pow(R, 2) + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
	temp += fp->mMassBulge / (2.0 * M_PI) * fp->aBulge / (r * pow(r + fp->aBulge, 3));

	return temp;
};

//what about black hole?!
double gslDensity(double z, void* p) { // function Units: SolarMassUnit*pc^-3

	struct gslDensityQagsParams* fp = (struct gslDensityQagsParams*)p;

	double z2b2 = pow(z, 2) + pow(fp->bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R2 = pow(fp->R, 2);

	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * R2 + (fp->aDisk + 3 * sz2b2) * pow(fp->aDisk + sz2b2, 2)) / (pow(R2 + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
	temp += fp->mMassBulge / (2 * M_PI) * fp->aBulge / (sqrt(R2 + pow(z, 2)) * pow(sqrt(R2 + pow(z, 2)) + fp->aBulge, 3));
	return temp;
};

//for MC
double gslDensityDisk(double x[], size_t dim, void* p) {
	struct gslDensityDiskParams* fp = (struct gslDensityDiskParams*)p;

	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double z2b2 = pow(x[2], 2) + pow(fp->bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R = gsl_hypot(x[0], x[1]);
	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * pow(R, 2) + (fp->aDisk + 3 * sz2b2) * pow(fp->aDisk + sz2b2, 2)) / (pow(pow(R, 2) + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

	return temp;
};

//for qags
double gslDensityDisk(double z, void* p) {
	struct gslDensityDiskQagsParams* fp = (struct gslDensityDiskQagsParams*)p;

	double z2b2 = pow(z, 2) + pow(fp->bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R2 = pow(fp->R, 2);
	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * R2 + (fp->aDisk + 3 * sz2b2) * pow(fp->aDisk + sz2b2, 2)) / (pow(R2 + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

	return temp;
};

double gslDensityBulge(double x[], size_t dim, void* p) {

	struct gslDensityBulgeParams* fp = (struct gslDensityBulgeParams*)p;

	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double r = gsl_hypot3(x[0], x[1], x[2]);
	double temp = fp->mMassBulge / (2 * M_PI * pow(fp->aBulge, 3)) * pow(fp->aBulge, 4) / (r * pow(r + fp->aBulge, 3));

	return temp;
};

double gslDensityBulge(double z, void* p) {

	struct gslDensityBulgeQagsParams* fp = (struct gslDensityBulgeQagsParams*)p;
	double R2 = pow(fp->R, 2);
	double temp = fp->mMassBulge / (2 * M_PI ) * fp->aBulge / (sqrt(R2 +pow(z,2)) * pow(sqrt(R2 + pow(z, 2)) + fp->aBulge, 3));

	return temp;
};

double Potential::closestToZero(double a, double b) {
	double temp = 0;
	if ((a < 0 && b> 0) || (a > 0 && b<0))
		temp = 0;
	if (a > 0 && b > 0) {
		if (a > b)
			return b;
		return a;
	}
	if (a < 0 && b < 0) {
		if (a > b)
			return a;
		return b;
	}
	return temp;
}

double Potential::circularVelocity(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y,2);
	double velocity = sqrt(G * (mMassBlackHole * R2 / pow(R2 + pow(z, 2), 1.5)
		+ mMassDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5)
		+ mMassBulge * R2 / (r * pow(aBulge + r, 2))
		+ 4 * M_PI * densityHalo * R2 * pow(rHalo, 3) * log(r / rHalo + 1) / pow(R2+pow(z,2), 1.5)
		- 4 * M_PI * densityHalo * R2 * pow(rHalo, 2) / (pow(r, 2) * (r / rHalo + 1))));
	return velocity;
}

double Potential::circularVelocityDisk(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt( G * mMassDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5));
	return velocity;
}

double Potential::circularVelocityBlackHole(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * mMassBlackHole * R2 / pow(R2 + pow(z, 2), 1.5));
	return velocity;
}

double Potential::circularVelocityBulge(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * mMassBulge * R2 / (r * pow(aBulge + r, 2)));
	return velocity;
}

double Potential::circularVelocityHalo(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(4 * M_PI * G * densityHalo * R2 * pow(rHalo, 3) * log(sqrt(R2 + pow(z, 2)) / rHalo + 1) / pow(R2 + pow(z, 2), 1.5)
		- 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 2) / ((R2+pow(z, 2)) * (sqrt(R2 + pow(z, 2)) / rHalo + 1)));
	return velocity;
}

double Potential::densityDisk(double R, double z){
	double z2b2 = pow(z, 2) + pow(bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	return pow(bDisk, 2)* mMassDisk / (4 * M_PI) * (aDisk * pow(R, 2) + (aDisk + 3 * sz2b2 * pow(aDisk + sz2b2, 2))) / (pow(pow(R, 2) + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
}

double Potential::densityDisk(double x, double y, double z){
	double R = sqrt(pow(x, 2) + pow(y, 2));
	return densityDisk(R,z);
}

double Potential::surfaceDensityDisk(double R){
	gsl_function F;
	F.function = &gslDensityDisk;
	gslDensityDiskQagsParams densityDiskQagsParams = { mMassDisk, aDisk, bDisk, R };
	F.params = &densityDiskQagsParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0*result;
}

Vec3D Potential::sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
	double smallestx = closestToZero(xMin, xMax);
	double smallesty = closestToZero(yMin, yMax);
	double smallestz = closestToZero(zMin, zMax);

	double k = densityDisk(smallestx, smallesty, smallestz);
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(xMin, xMax);
	std::uniform_real_distribution<> disy(yMin, yMax);
	std::uniform_real_distribution<> disz(zMin, zMax);
	std::uniform_real_distribution<> disaccept(0, k);

	while (true) {
		double x = disx(gen);
		double y = disy(gen);
		double z = disz(gen);

		double accept = disaccept(gen);
		double temp = densityDisk(x, y, z);

		if (accept < temp)
			return Vec3D(x, y, z);
	}
	return Vec3D();
}


double Potential::densityBulge(double r){
	return mMassBulge/(2*M_PI*pow(aBulge,3))*pow(aBulge,4)/(r*pow(r+aBulge,3));
}

double Potential::densityBulge(double x, double y, double z){
	double r = sqrt(x * x + y * y + z * z);
	return densityBulge(r);
}

double Potential::surfaceDensityBulge(double R){
	if (R == 0) {
		std::cout << "R=0 not allowed." << std::endl;
		return 0;
	}
	//1/(2*M_PI)*1/pow(1-pow(R,2))^2*((2+pow(R,2)*gsl_arcsec gsl_complex_arcsec_real(R)/sqrt(pow(R,2)-1))
	gsl_function F;
	F.function = &gslDensityBulge;
	gslDensityBulgeQagsParams densityBulgeQagsParams = { mMassBulge,aBulge,R };
	F.params = &densityBulgeQagsParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F,0,0,1e-7, 1000,iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0 * result;
}

//todo: change to vector<Star> as soon as present day mass function is sorted out
Vec3D Potential::sampleBuldge(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
	double k = 0;
	double smallestx = closestToZero(xMin, xMax);
	double smallesty = closestToZero(yMin, yMax);
	double smallestz = closestToZero(zMin, zMax);
	if (sqrt(smallestx * smallestx + smallesty * smallesty + smallestz * smallestz < 0.1)) {
		//std::cout << "Warning: Only samling Stars for r>100pc" << std::endl;
		k = densityBulge(0.1);
	}
	else
		k = densityBulge(smallestx, smallesty, smallestz);
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> disx(xMin, xMax);
	std::uniform_real_distribution<> disy(yMin, yMax);
	std::uniform_real_distribution<> disz(zMin, zMax);
	std::uniform_real_distribution<> disaccept(0, k);

	while (true) {
		double x = disx(gen);
		double y = disy(gen);
		double z = disz(gen);
		if (sqrt(x * x + y * y + z * z) < 0.1)
			continue;

		double accept = disaccept(gen);
		double temp = densityBulge(x, y, z);

		if (accept < temp)
			return Vec3D(x, y, z);
	}
	return Vec3D();
}

double Potential::frequencyDistribution(Vec3D position, Vec3D volumeElement) {

	struct gslDensityParams frequencyDistributionParams = { mMassBulge, aBulge, mMassDisk, aDisk, bDisk };

	gsl_monte_function F;

	F.f = &gslDensity;
	F.dim = 3;
	F.params = &frequencyDistributionParams;

	gsl_rng_env_setup();

	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double xl[3] = { position.x, position.y, position.z };
	double xu[3] = { position.x+volumeElement.x, position.y+volumeElement.y, position.z+volumeElement.z };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double res, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,
		&res, &err);
	gsl_monte_plain_free(s);

	return res;
}

double Potential::massDisk(Vec3D position, Vec3D volumeElement){

	struct gslDensityDiskParams densityDiskParams = { mMassDisk, aDisk, bDisk };

	gsl_monte_function F;

	F.f = &gslDensityDisk;
	F.dim = 3;
	F.params = &densityDiskParams;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double xl[3] = { position.x, position.y, position.z };
	double xu[3] = { position.x + volumeElement.x, position.y + volumeElement.y, position.z + volumeElement.z };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double mass, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,
		&mass, &err);
	gsl_monte_plain_free(s);

	return mass;
}

double Potential::massBulge(Vec3D position, Vec3D volumeElement) {

	struct gslDensityBulgeParams densityBulgeParams = { mMassBulge, aBulge };

	gsl_monte_function F;

	F.f = &gslDensityBulge;
	F.dim = 3;
	F.params = &densityBulgeParams;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double xl[3] = { position.x, position.y, position.z };
	double xu[3] = { position.x + volumeElement.x, position.y + volumeElement.y, position.z + volumeElement.z };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double mass, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,
		&mass, &err);
	gsl_monte_plain_free(s);

	return mass;
}

double Potential::angularVelocity(double R){
	double R2 = gsl_pow_2(R);
	double haloTemp = 4.0 * M_PI * densityHalo * gsl_pow_3(rHalo);
	return kmInpc*sqrt(G/R*(-mMassBulge/ gsl_pow_2(aBulge+R)+mMassDisk*R/pow(gsl_pow_2(aDisk+bDisk)+R2,1.5)
		+mMassBlackHole/R2+ haloTemp/(R2+R*rHalo)+ haloTemp*log((R+rHalo)/rHalo)/R2));
}

double Potential::surfaceDensity(double R){
	if (R == 0) {
		std::cout << "R=0 not allowed." << std::endl;
		return 0;
	}
	
	gsl_function F;
	F.function = &gslDensity;
	gslDensityQagsParams densityQagsParams = { mMassBulge, aBulge, mMassDisk, aDisk, bDisk, R };
	F.params = &densityQagsParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0 * result;
}

double Potential::epicyclicFrequency(double R, double z){
	double r = gsl_hypot(R, z);
	double r2 = gsl_pow_2(r);
	double r3 = gsl_pow_3(r);
	double r4 = gsl_pow_2(r2);
	double r5 = gsl_pow_5(r);
	double R2 = gsl_pow_2(R);
	//double R3 = gsl_pow_3(R);
	//double R4 = gsl_pow_2(R2);
	//double R6 = gsl_pow_2(R3);
	//double haloTemp = 4.0 * M_PI * densityHalo * gsl_pow_3(rHalo);
	//double temp = mMassBlackHole / R3 + 2.0 * mMassBulge / gsl_pow_3(aBulge + R) - 3.0 * mMassBulge / (R * gsl_pow_2(aBulge + R));
	//temp = temp - 3.0 * mMassDisk * R2 / pow(gsl_pow_2(aDisk) + 2.0 * aDisk * bDisk + gsl_pow_2(bDisk) + R2, 2.5) + 4.0 * mMassDisk / pow(gsl_pow_2(aDisk) + 2.0 * aDisk * bDisk + gsl_pow_2(bDisk) + R2, 1.5);
	//temp = temp + haloTemp / (R * gsl_pow_2(R + rHalo)) - haloTemp / (R2 * (R + rHalo)) + haloTemp * log((R + rHalo) / rHalo) / R3;
	double temp = -3.0 * mMassBlackHole * R2 / r5 + 4.0 * mMassBlackHole / gsl_pow_3(r);
	temp += 4.0 * mMassDisk / pow(gsl_pow_2(aDisk) + R2 + 2.0 * aDisk * sqrt(gsl_pow_2(bDisk) + gsl_pow_2(z)) + gsl_pow_2(bDisk) + gsl_pow_2(z), 1.5) - 3.0* mMassDisk * R2 / pow(gsl_pow_2(aDisk) + R2 + 2.0* aDisk*sqrt(gsl_pow_2(bDisk) + gsl_pow_2(z)) + gsl_pow_2(bDisk) + gsl_pow_2(z),2.5);
	temp += 2.0 * mMassBulge * R2 / (r2 * gsl_pow_3(aBulge + r)) + (mMassBulge * R2) / (gsl_pow_3(r) * gsl_pow_2(aBulge + r)) - (4.0 * mMassBulge) / gsl_pow_2(r * (aBulge + r));
	temp += 12.5664 * densityHalo * R2 * gsl_pow_3(rHalo) / (gsl_pow_3(r) * gsl_pow_2(rHalo + r)) + 37.6991 * densityHalo * R2 * gsl_pow_3(rHalo) / (r4 * (rHalo + r)) - 50.2655 * densityHalo * gsl_pow_3(rHalo) / (r2 * (rHalo + r)) + densityHalo * gsl_pow_3(rHalo) * (50.2655 / r3 - 37.6991 * R2 / r5) * log((rHalo + r) / rHalo);

	if (temp < 0) {
		std::cout << "Warning: Imaginary epicyclicFrequency. Setting 1e-10 to avoid crash" << std::endl;
		temp = 1e-10;
	}
	temp = sqrt(G * temp)*kmInpc;
	return temp;
}

double Potential::radialVelocityDispersion(double R, double z){
	//3.36 * G * surfaceDensity(R) / epicyclicFrequency(R)*kmInpc;
	//assumed to be normalized ...
	double temp = exp(-gsl_pow_2(R) + 2 * gsl_pow_2(velocityDispersionScaleLength) / aDisk) * 3.36 * G * surfaceDensity(R) / epicyclicFrequency(R,z) * kmInpc;
	if (temp < 0) {
		std::cout << "Warning: Negative radial velocity dispersion. Setting 1 to avoid crash" << std::endl;
		return 1;
	}
	return temp;
}

double Potential::verticalVelocityDispersion(double R){
	return sqrt(M_PI* G* surfaceDensity(R)* bDisk);
}

double Potential::azimuthalVelocityDispersion(double R, double z){
	return radialVelocityDispersion(R,z)* gsl_pow_2(epicyclicFrequency(R,z)) / (4.0 * gsl_pow_2(angularVelocity(R)));
}

double Potential::azimuthalStreamingVelocity(Vec3D position){
	double R = gsl_hypot(position.x,position.y);
	return sqrt(gsl_pow_2(radialVelocityDispersion(R,position.z)) * (1.0 - gsl_pow_2(epicyclicFrequency(R, position.z)) / (4.0 * gsl_pow_2(angularVelocity(R))) - 2.0 * R / aDisk) + gsl_pow_2(circularVelocity(&position)));
}

void Potential::applyForce(Star* star){
	double r = star->position.length();
	double r2 = gsl_pow_2(r);
	double r3 = gsl_pow_3(r);
	double tempHalo = 4.0 * M_PI * densityHalo;
	//BH, Halo & Bulge
	double temp = -G * (mMassBlackHole / r3 + mMassBulge / (r * gsl_pow_2(aBulge + r)) - tempHalo * gsl_pow_2(rHalo) / (r2 * (1 + r / rHalo)) + tempHalo * gsl_pow_3(rHalo) * log(1 + r / rHalo)/ r3)*kmInpc;
	star->acceleration.x += temp * star->position.x;
	star->acceleration.y += temp * star->position.y;
	star->acceleration.z += temp * star->position.z;
	//Disk
	double x2y2 = gsl_pow_2(star->position.x) + gsl_pow_2(star->position.y);
	double sqrtb2z2 = gsl_hypot(bDisk, star->position.z);
	temp = -G * mMassDisk / pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5)*kmInpc;
	star->acceleration.x += temp * star->position.x;
	star->acceleration.y += temp * star->position.y;
	star->acceleration.z += -G * mMassDisk * star->position.z * (aDisk + sqrtb2z2) / (sqrtb2z2 * pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5)) * kmInpc;
}
