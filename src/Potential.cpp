#include "..\include\Potential.h"

struct gslDensityDiskParams { double mMassDisk; double aDisk; double bDisk; };
struct gslDensityBulgeParams { double mMassBulge; double aBulge; };
struct gslDensityParams { double mMassBulge; double aBulge; double mMassDisk; double aDisk; double bDisk; };

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

Potential::Potential(Vec3D position){
	this->position = position;
}


double Potential::circularVelocity(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y,2);
	double velocity = sqrt(G * massBlackHole * R2 / pow(R2 + pow(z, 2), 1.5)
		+ G * mMassDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5)
		+ G * mMassBulge * R2 / (r * pow(aBulge + r, 2))
		+ 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 3) * log(r / rHalo + 1) / pow(R2+pow(z,2), 1.5)
		- 4 * M_PI * G * densityHalo * R2 * pow(rHalo, 2) / (pow(r, 2) * (r / rHalo + 1)));
	return velocity;
}

double Potential::circularVelocityDisk(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt( G * mMassDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5));
	return velocity;
}

double Potential::circularVelocityBlackHole(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * massBlackHole * R2 / pow(R2 + pow(z, 2), 1.5));
	return velocity;
}

double Potential::circularVelocityBulge(Vec3D* position){
	Vec3D distance = *position - this->position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(G * mMassBulge * R2 / (r * pow(aBulge + r, 2)));
	return velocity;
}

double Potential::circularVelocityHalo(Vec3D* position){
	Vec3D distance = *position - this->position;
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


double gslDensity(double x[], size_t dim, void* p) {

	struct gslDensityParams* fp = (struct gslDensityParams*)p;

	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double z2b2 = pow(x[2], 2) + pow(fp->bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R = gsl_hypot(x[0],x[1]);
	double r = gsl_hypot3(x[0], x[1], x[2]);
	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * pow(R, 2) + (fp->aDisk + 3 * sz2b2 * pow(fp->aDisk + sz2b2, 2))) / (pow(pow(R, 2) + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
	temp += fp->mMassBulge / (2 * M_PI * pow(fp->aBulge, 3)) * pow(fp->aBulge, 4) / (r * pow(r + fp->aBulge, 3));

	return temp;
};

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
	double temp = gsl_pow_2(fp->bDisk) * fp->mMassDisk / (4 * M_PI) * (fp->aDisk * pow(R, 2) + (fp->aDisk + 3 * sz2b2 * pow(fp->aDisk + sz2b2, 2))) / (pow(pow(R, 2) + pow(fp->aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

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