#include "..\include\MWPotential.h"

const double MWPotential::mMassBlackHole = 4e6; // SolarMassUnit
const double MWPotential::mMassDisk = 6.51e10; // SolarMassUnit
const double MWPotential::aDisk = 4.4e3; // pc
const double MWPotential::bDisk = 0.267e3; // pc
const double MWPotential::mMassBulge = 1.8e10; // SolarMassUnit
const double MWPotential::mMassSmallBulge = 0; // SolarMassUnit

Hernquist MWPotential::bulgePotential = Hernquist(mMassBulge, aBulge);

const double MWPotential::aBulge = 0.9e3;//0.35e3; // pc
const double MWPotential::aSmallBulge = 0.35e3; // pc
const double MWPotential::characteristicVelocityBulge = 444.4;  // km/s
const double MWPotential::rHalo = 17e3; // pc
const double MWPotential::mMassHalo = 5E11; // SolarMassUnit
//const double MWPotential::kmInpc = 3.086e-13;
const double MWPotential::velocityDispersionScaleLength = aDisk / 4.0;

const std::string MWPotential::lookupTableLocation = "src/LookupTables/";
const std::string MWPotential::velocityDistributionBulgeTableFilename = "velocityDistributionBulgeTable.dat";

struct gslRParam { double R; };
struct gslDensityConeParams { Matrix* transformation; double distance; double r; double x; double y; };

struct gslSphericalAveragedDiskParams { double mMassDisk; double aDisk; double bDisk; double r; };

double MWPotential::gslDensity(double x[], size_t dim, void* p) {

	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double R = gsl_hypot(x[0], x[1]);
	double z = x[2];

	return densityDisk(R, z) + bulgePotential.density(R, z);
};

double MWPotential::gslDensity(double z, void* p) { // function Units: SolarMassUnit*pc^-3

	struct gslRParam* fp = (struct gslRParam*)p;

	double temp = MWPotential::densityDisk(fp->R,z) + bulgePotential.density(fp->R, z);
	return temp;
};

//for Monte Carlo Integration
double MWPotential::gslDensityDisk(double x[], size_t dim, void* p) {
	if (dim != 3)
	{
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	double z2b2 = pow(x[2], 2) + pow(bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R = gsl_hypot(x[0], x[1]);
	double temp = gsl_pow_2(bDisk) * mMassDisk / (4 * M_PI) * (aDisk * pow(R, 2) + (aDisk + 3 * sz2b2) * pow(aDisk + sz2b2, 2)) / (pow(pow(R, 2) + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

	return temp;
};

//for qags
double MWPotential::gslDensityDisk(double z, void* p) {
	struct gslRParam* fp = (struct gslRParam*)p;

	double z2b2 = pow(z, 2) + pow(bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	double R2 = pow(fp->R, 2);
	double temp = gsl_pow_2(bDisk) * mMassDisk / (4 * M_PI) * (aDisk * R2 + (aDisk + 3 * sz2b2) * pow(aDisk + sz2b2, 2)) / (pow(R2 + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

	return temp;
}

double MWPotential::gslDensityDiskx(double x, void* p){
	gsl_function F;
	F.function = &gslDensityDisky;
	struct gslDensityConeParams* fp = (struct gslDensityConeParams*)p;
	fp->x = x;
	F.params = fp;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;
	double boundary = sqrt(fp->r * fp->r - x * x);
	gsl_integration_qag(&F, -boundary, boundary, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return result;
}

double MWPotential::gslDensityDisky(double y, void* p){
	gsl_function F;
	F.function = &gslDensityDiskz;
	struct gslDensityConeParams* fp = (struct gslDensityConeParams*)p;
	fp->y = y;
	F.params = fp;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_integration_qag(&F, fp->distance/fp->r*sqrt(fp->x*fp->x+y*y), fp->distance, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return result;
}

double MWPotential::gslDensityDiskz(double z, void* p){
	struct gslDensityConeParams* fp = (struct gslDensityConeParams*)p;
	Vec3D location = Vec3D(fp->x, fp->y, z);
	location = *(fp->transformation) * location;

	double z2b2 = location.z*location.z + bDisk*bDisk;
	double sz2b2 = sqrt(z2b2);
	double R = gsl_hypot(location.x, location.y);
	double temp = gsl_pow_2(bDisk) * mMassDisk / (4 * M_PI) * (aDisk * pow(R, 2) + (aDisk + 3 * sz2b2) * pow(aDisk + sz2b2, 2)) / (pow(pow(R, 2) + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));

	return temp;
}

double MWPotential::gslVelocityBulge(double r, void* p){
	double r2 = gsl_pow_2(r);
	//struct gslVelocityBulgeParams* fp = (struct gslVelocityBulgeParams*)p;
	//double temp = 1/(r * gsl_pow_3(fp->aBulge + r)) * (
	double temp = mMassBlackHole / r2 + bulgePotential.potentialdr(r);
		//- mMassHalo / (r2 + r * rHalo) + mMassHalo * log((r + rHalo) / rHalo) / r2;
	if (r < aBulge)
		temp += MWPotential::sphericalAveragedDisc(r);
	temp = temp * MWPotential::bulgePotential.density(r);
	//double temp = 1 / (r * gsl_pow_3(fp->aBulge + r)) * (fp->mMassBlackHole / r2 + fp->mMassBulge / gsl_pow_2(fp->aBulge + r) + fp->mMassDisk * r / pow(gsl_pow_2(fp->aDisk + fp->bDisk) + r2, 1.5) - haloTemp / (r2 + r * fp->rHalo) + haloTemp * log((r + fp->rHalo) / fp->rHalo) / r2);
	//double temp = 1 / (r * gsl_pow_3(fp->aBulge + r)) * (fp->mMassBlackHole / r2 + fp->mMassBulge / gsl_pow_2(fp->aBulge + r) + fp->mMassDisk * r *(1/3*fp->aDisk+0.7698*sqrt(3*fp->bDisk+r2))/ (sqrt(gsl_pow_2(fp->bDisk)+r2/3)*pow(r2+gsl_pow_2(fp->aDisk+gsl_pow_2(fp->aDisk+sqrt(gsl_pow_2(fp->bDisk)+r2/3))),1.5)) - haloTemp / (r2 + r * fp->rHalo) + haloTemp * log((r + fp->rHalo) / fp->rHalo) / r2);
	//double temp = 1 / (r * gsl_pow_3(fp->aBulge + r)) * (fp->mMassBulge / gsl_pow_2(fp->aBulge + r));
	return temp;
}

double gslSphericalAveragedDisk(double theta, void* p) {
	struct gslSphericalAveragedDiskParams* fp = (struct gslSphericalAveragedDiskParams*)p;
	double temp = gsl_hypot(fp->bDisk, fp->r * cos(theta));
	//return 0.5 * fp->mMassDisk * (2. * fp->r * gsl_pow_2(cos(theta)) * (fp->aDisk + temp) / temp + 2 * fp->r * gsl_pow_2(sin(theta))) / pow(gsl_pow_2(fp->aDisk + temp) + gsl_pow_2(fp->r * sin(theta)), 1.5);
	double r2 = gsl_pow_2(fp->r);
	return fp->mMassDisk / pow(r2 * gsl_pow_2(sin(theta)) + gsl_pow_2(fp->aDisk + gsl_hypot(fp->r * cos(theta), fp->bDisk)), 1.5);
}

double MWPotential::closestToZero(double a, double b) {
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

MWPotential::MWPotential(): velocityDistributionBulgeTable(LookupTable(velocityDistributionBulgeTableFilename)) {
	if (velocityDistributionBulgeTable.isEmpty()) {
		std::cout << "No lookup table found for the bulge velocity distribution. Generating new table ... " << std::endl;
		generateVelocityDistributionBulgeLookupTable(25000);
		velocityDistributionBulgeTable = LookupTable(velocityDistributionBulgeTableFilename);
	}
	this->bulgePotential = Hernquist(mMassBulge, aBulge);
}

double MWPotential::sphericalAveragedDisc(double r){
	//gslSphericalAveragedDisk(double theta, void* p) {
	//struct gslSphericalAveragedDiskParams* fp = (struct gslSphericalAveragedDiskParams*)p;

	gsl_function F;
	F.function = &gslSphericalAveragedDisk;
	gslSphericalAveragedDiskParams sphericalAveragedDiskParams = {mMassDisk, aDisk, bDisk, r };
	F.params = &sphericalAveragedDiskParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qag(&F, 0, 2*M_PI, 1e-3, 1e-3, 1000,1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 1/(2*M_PI) * result;
}

double MWPotential::circularVelocity(const Vec3D* position){
	double r = position->length();
	double r2 = gsl_pow_2(r);
	double r3 = r2 * r;
	double z = position->z;
	double R2 = pow(position->x, 2) + pow(position->y,2);
	double velocity = sqrt(Constants::G * R2*mMassBlackHole/ r3
		+ Constants::G * R2 * mMassDisk/ pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5)
		+ gsl_pow_2(bulgePotential.circularVelocity(position))
		+ Constants::G * R2 * mMassHalo * log(r / rHalo + 1) / r3
		- Constants::G * R2 * mMassHalo / (r2 * (rHalo  + r)));
	return velocity;
}

double MWPotential::circularVelocityDisk(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(Constants::G * mMassDisk * R2 / pow(pow(aDisk + sqrt(pow(bDisk, 2) + pow(z, 2)), 2) + R2, 1.5));
	return velocity;
}

double MWPotential::circularVelocityBlackHole(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(Constants::G * mMassBlackHole * R2 / pow(R2 + pow(z, 2), 1.5));
	return velocity;
}

double MWPotential::circularVelocityHalo(Vec3D* position){
	Vec3D distance = *position;
	double r = distance.length();
	double r2 = gsl_pow_2(r);
	double z = distance.z;
	double R2 = pow(distance.x, 2) + pow(distance.y, 2);
	double velocity = sqrt(Constants::G * R2 * mMassHalo * ( log(r / rHalo + 1) / (r*r2) - 1 / (r2 * (r + rHalo))));
	return velocity;
}

double MWPotential::escapeVelocity(Vec3D* position){
	double r = position->length();
	double R = gsl_hypot(position->x, position->y);
	return sqrt(2* Constants::G *mMassBlackHole/r+ 2 * Constants::G * mMassDisk/sqrt(gsl_pow_2(R)+gsl_pow_2(aDisk+gsl_hypot(position->z,bDisk)))+ gsl_pow_2(bulgePotential.escapeVelocity(r))+ 2 * Constants::G * mMassHalo*log(1+r/rHalo)/r);
}

double MWPotential::densityDisk(double R, double z){
	double R2 = gsl_pow_2(R);
	double z2b2 = pow(z, 2) + pow(bDisk, 2);
	double sz2b2 = sqrt(z2b2);
	return pow(bDisk, 2)* mMassDisk / (4 * M_PI) * (aDisk * R2 + (aDisk + 3 * sz2b2) * pow(aDisk + sz2b2, 2)) / (pow(R2 + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
}

double MWPotential::densityDisk(double x, double y, double z){
	double R = sqrt(pow(x, 2) + pow(y, 2));
	return densityDisk(R,z);
}

double MWPotential::denistyWang(double x, double y, double z){
	double density0 = 2.1242;
	double q = 0.6;
	double z0 = 0.4e3;
	double x0 = 1.49e3;
	double y0 = 0.58e3;
	double r1 = pow(gsl_pow_2(gsl_pow_2(x / x0) + gsl_pow_2(y / y0)) + gsl_pow_4(z / z0), 0.25);
	double r2 = sqrt((gsl_pow_2(q) * (gsl_pow_2(x) + gsl_pow_2(y)) + gsl_pow_2(z)) / gsl_pow_2(z0));
	return density0 * (exp(-gsl_pow_2(r1) / 2) + pow(r2, -1.85) * exp(-r2));
}

double MWPotential::surfaceDensityDisk(double R){
	gsl_function F;
	F.function = &gslDensityDisk;
	gslRParam densityDiskQagsParams = { R };
	F.params = &densityDiskQagsParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0*result;
}

double MWPotential::massDisk(Vec3D position, Vec3D volumeElement){
	gsl_monte_function F;

	F.f = &gslDensityDisk;
	F.dim = 3;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double x1_low = std::min(position.x, position.x + volumeElement.x);
	double x1_high = std::max(position.x, position.x + volumeElement.x);
	double x2_low = std::min(position.y, position.y + volumeElement.y);
	double x2_high = std::max(position.y, position.y + volumeElement.y);
	double x3_low = std::min(position.z, position.z + volumeElement.z);
	double x3_high = std::max(position.z, position.z + volumeElement.z);
	double xl[3] = { x1_low, x2_low, x3_low };
	double xu[3] = { x1_high, x2_high, x3_high };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double mass, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s, &mass, &err);
	gsl_monte_plain_free(s);

	return mass;
}

double MWPotential::massDisk(Matrix* transformation, double distance, double r) {
	gsl_function F;
	F.function = &gslDensityDiskx;
	gslDensityConeParams densityDiskConeParam{ transformation, distance, r,0,0 };
	F.params = &densityDiskConeParam;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qag(&F, -r, r, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return result;
}

double MWPotential::angularVelocity(double R){
	double R2 = gsl_pow_2(R);
	return Constants::kmInpc*sqrt(Constants::G /R*(bulgePotential.potentialdr(R)+mMassDisk*R/pow(gsl_pow_2(aDisk+bDisk)+R2,1.5)
		+mMassBlackHole/R2+ mMassHalo/(R2+R*rHalo)+ mMassHalo*log((R+rHalo)/rHalo)/R2));
}

double MWPotential::surfaceDensity(double R){
	if (R == 0) {
		std::cout << "R=0 not allowed." << std::endl;
		return 0;
	}
	
	gsl_function F;
	F.function = &MWPotential::gslDensity;
	gslRParam densityQagsParams = {R};
	F.params = &densityQagsParams;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0 * result;
}

double MWPotential::epicyclicFrequency(double R, double z){
	double r = gsl_hypot(R, z);
	double r2 = gsl_pow_2(r);
	double r3 = gsl_pow_3(r);
	double r4 = gsl_pow_2(r2);
	double r5 = gsl_pow_5(r);
	double R2 = gsl_pow_2(R);

	double temp = -3.0 * mMassBlackHole * R2 / r5 + 4.0 * mMassBlackHole / r3;
	temp += 4.0 * mMassDisk / pow(gsl_pow_2(aDisk) + R2 + 2.0 * aDisk * sqrt(gsl_pow_2(bDisk) + gsl_pow_2(z)) + gsl_pow_2(bDisk) + gsl_pow_2(z), 1.5) - 3.0* mMassDisk * R2 / pow(gsl_pow_2(aDisk) + R2 + 2.0* aDisk*sqrt(gsl_pow_2(bDisk) + gsl_pow_2(z)) + gsl_pow_2(bDisk) + gsl_pow_2(z),2.5);
	temp += bulgePotential.potentialdr2(R)+3./R*bulgePotential.potentialdr(R);
	temp += mMassHalo*R2 / (r3 * gsl_pow_2(rHalo + r)) + 3 * R2* mMassHalo / (r4 * (rHalo + r)) - 4*mMassHalo / (r2 * (rHalo + r)) + mMassHalo * (4 / r3 - 3 * R2 / r5) * log((rHalo + r) / rHalo);

	if (temp < 0) {
		std::cout << "Warning: Imaginary epicyclicFrequency. Setting 1e-10 to avoid crash" << std::endl;
		temp = 1e-10;
	}
	temp = sqrt(Constants::G * temp)* Constants::kmInpc;
	return temp;
}

double MWPotential::radialVelocityDispersionDisk(double R, double z){
	//3.36 * G * surfaceDensity(R) / epicyclicFrequency(R)*kmInpc;
	//assumed to be normalized ...
	double temp = exp(-gsl_pow_2(R) + 2 * gsl_pow_2(velocityDispersionScaleLength) / aDisk) * 3.36 * Constants::G * surfaceDensity(R) / epicyclicFrequency(R,z) * Constants::kmInpc;
	if (temp < 0) {
		std::cout << "Warning: Negative radial velocity dispersion. Setting 1 to avoid crash" << std::endl;
		return 1;
	}
	return temp;
}

double MWPotential::verticalVelocityDispersion(double R){
	return sqrt(M_PI* Constants::G * surfaceDensity(R)* bDisk);
}

double MWPotential::azimuthalVelocityDispersion(double R, double z){
	return radialVelocityDispersionDisk(R,z)* gsl_pow_2(epicyclicFrequency(R,z)) / (4.0 * gsl_pow_2(angularVelocity(R)));
}

double MWPotential::azimuthalStreamingVelocity(Vec3D position){
	double R = gsl_hypot(position.x,position.y);
	return sqrt(gsl_pow_2(radialVelocityDispersionDisk(R,position.z)) * (1.0 - gsl_pow_2(epicyclicFrequency(R, position.z)) / (4.0 * gsl_pow_2(angularVelocity(R))) - 2.0 * R / aDisk) + gsl_pow_2(circularVelocity(&position)));
}

double MWPotential::velocityDispersionBulge(double r){
	if (r == 0 && debug) {
		std::cout << "\nr=0 not allowed." << std::endl;
		return 0;
	}
	if (r == 0) 
		return 0;

	gsl_function F;
	F.function = &gslVelocityBulge;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, r, 0, 1e-5, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return sqrt(Constants::G /MWPotential::bulgePotential.density(r) * result);
}

double MWPotential::velocityDistributionBulgeTableValue(double r){
	return velocityDistributionBulgeTable.get(r);
}

void MWPotential::generateVelocityDistributionBulgeLookupTable(double rMax){
	LookupTable lookupTable = LookupTable(velocityDistributionBulgeTableFilename);
	std::string header = "mMassBlackHole " + std::to_string(mMassBlackHole) +
		", aBulge " + std::to_string(aBulge) +", mMassDisk " + std::to_string(mMassDisk) +", mMassBulge " + std::to_string(mMassBulge) +
		", aDisk " + std::to_string(aDisk) +", bDisk " + std::to_string(bDisk) +", rHalo " + std::to_string(rHalo) +", mMassHalo " + std::to_string(mMassHalo);
	std::vector<double> rVec;
	std::vector<double> velDist;
	ProgressBar progressBar = ProgressBar(0, rMax);
	for (double r = 1; r <= rMax; ++r) {
		progressBar.Update(r);
		progressBar.Print();
		rVec.push_back(r);
		velDist.push_back(velocityDispersionBulge(r));
	}
	lookupTable.setMap(rVec, velDist);
	lookupTable.makeFile(velocityDistributionBulgeTableFilename, header);
	return;
}

void MWPotential::applyForce(Star* star){
	double r = star->position.length();
	double r2 = r*r;
	double r3 = r2*r;
	//BH, Halo & Bulge
	double temp = -Constants::G * (mMassBlackHole / r3 + bulgePotential.forceTemp(r) - mMassHalo / (rHalo*r2 * (1 + r / rHalo)) + mMassHalo * log(1 + r / rHalo)/ r3)* Constants::kmInpc;
	star->acceleration.x += temp * star->position.x;
	star->acceleration.y += temp * star->position.y;
	star->acceleration.z += temp * star->position.z;
	//Disk
	double x2y2 = gsl_pow_2(star->position.x) + gsl_pow_2(star->position.y);
	double sqrtb2z2 = gsl_hypot(bDisk, star->position.z);
	temp = -Constants::G * mMassDisk / pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5)* Constants::kmInpc;
	star->acceleration.x += temp * star->position.x;
	star->acceleration.y += temp * star->position.y;
	star->acceleration.z += -Constants::G * mMassDisk * star->position.z * (aDisk + sqrtb2z2) / (sqrtb2z2 * pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5)) * Constants::kmInpc;
}

void MWPotential::applyForce(Vec3D position, Vec3D& acceleration){
	double r = position.length();
	double r2 = r * r;
	double r3 = r2 * r;
	//BH, Halo & Bulge
	double temp = -Constants::G * (mMassBlackHole / r3 + bulgePotential.forceTemp(r) - mMassHalo / (rHalo * r2 * (1 + r / rHalo)) + mMassHalo * log(1 + r / rHalo) / r3) * Constants::kmInpc;
	acceleration.x += temp * position.x;
	acceleration.y += temp * position.y;
	acceleration.z += temp * position.z;
	//Disk
	double x2y2 = gsl_pow_2(position.x) + gsl_pow_2(position.y);
	double sqrtb2z2 = gsl_hypot(bDisk, position.z);
	temp = -Constants::G * mMassDisk / pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5) * Constants::kmInpc;
	acceleration.x += temp * position.x;
	acceleration.y += temp * position.y;
	acceleration.z += -Constants::G * mMassDisk * position.z * (aDisk + sqrtb2z2) / (sqrtb2z2 * pow(x2y2 + gsl_pow_2(aDisk + sqrtb2z2), 1.5)) * Constants::kmInpc;
}












//struct gslVelocityDispersionBulgeParams { double mMassBlackHole; double mMassBulge; double aBulge; double mMassDisk; double aDisk; double bDisk; double rHalo; double mMassHalo; double G; double R; };

//double MWPotential::radialVelocityDispersionBulge(double R, double z){
//	double r = gsl_hypot(R, z);
//	if (r == 0) {
//		std::cout << "r=0 not allowed." << std::endl;
//		return 0;
//	}
//
//	gsl_function F;
//	F.function = &gslVelocityDispersionBulge;
//	gslVelocityDispersionBulgeParams velocityDispersionBulgeParams = { mMassBlackHole, mMassBulge, aBulge, mMassDisk, aDisk, bDisk, rHalo, mMassHalo, G, R };
//	F.params = &velocityDispersionBulgeParams;
//	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
//	double result, error;
//
//	gsl_integration_qagiu(&F, abs(z), 0, 1e-7, 1000, iw, &result, &error);
//	gsl_integration_workspace_free(iw);
//	//return sqrt(r * gsl_pow_3(aBulge + r)* result);
//	return sqrt(1/(aBulge*mMassBulge/(2*M_PI*r*gsl_pow_3(aBulge+r))+densityHalo*gsl_pow_3(rHalo)/(r*gsl_pow_2(rHalo+r))+gsl_pow_2(bDisk)*mMassDisk*(aDisk*gsl_pow_2(R)+gsl_pow_2(aDisk+gsl_hypot(bDisk,z))*(aDisk + 3*gsl_hypot(bDisk, z)))/(4*M_PI*pow(gsl_pow_2(z)+gsl_pow_2(bDisk),1.5)*pow(gsl_pow_2(R)+ gsl_pow_2(aDisk + gsl_hypot(bDisk, z)),2.5)))  * result);
//}

//double MWPotential::infiniteDistributionFunctionBulge(double q){
//	double q2 = gsl_pow_2(q);
//	return 1/(pow(2,3.5)*gsl_pow_3(M_PI))*1/pow(1-q2,2.5)*(3*asin(q)+q*sqrt(1-q2)*(1-2*q2)*(8*gsl_pow_2(q2)-8*q2-3));
//}
//
//double MWPotential::distributionFunctionBulge(double e){
//	double q = sqrt(e) / characteristicVelocityBulge;
//	double qb = 0.787;
//	if (q < qb || q>=1)
//		return 0;
//	//1/(gsl_pow_2(aBulge2)*characteristicVelocityBulge)*
//	return (infiniteDistributionFunctionBulge(q)-infiniteDistributionFunctionBulge(qb));
//}
//
//double MWPotential::particleEnergy(Star* star){
//	return particleEnergy(star->position,star->velocity);
//}
//
//double MWPotential::particleEnergy(Vec3D& position, Vec3D& velocity){
//	return particleEnergy(position,velocity.length());
//}
//
//double MWPotential::particleEnergy(Vec3D& position, double velocity){
//	return 0.5 * gsl_pow_2(velocity) + potentialEnergy(position);
//}


//double MWPotential::potentialEnergy(Vec3D& position){
//	return MWPotential::potentialEnergy(gsl_hypot(position.x,position.y), position.z);
//}

//double MWPotential::potentialEnergy(double R, double z){
//	double r = gsl_hypot(R, z);
//	return -G*(mMassBulge2 / (r + aBulge2));
//	//return -G * (mMassBlackHole / r + mMassDisk / sqrt(gsl_pow_2(R) + gsl_pow_2(aDisk + sqrt(gsl_pow_2(z) + gsl_pow_2(bDisk))))+mMassBulge2/(r+aBulge2)+4*M_PI*gsl_pow_3(rHalo)*log(1+r/rHalo)/r);
//}

//double MWPotential::frequencyDistribution(Vec3D position, Vec3D volumeElement) {
//	gsl_monte_function F;
//
//	F.f = &gslDensity;
//	F.dim = 3;
//
//	gsl_rng_env_setup();
//
//	const gsl_rng_type *T = gsl_rng_default;
//	gsl_rng* r = gsl_rng_alloc(T);
//
//	double xl[3] = { position.x, position.y, position.z };
//	double xu[3] = { position.x+volumeElement.x, position.y+volumeElement.y, position.z+volumeElement.z };
//	size_t calls = 500000;
//
//	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
//	double res, err;
//	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,&res, &err);
//	gsl_monte_plain_free(s);
//	gsl_rng_free(r);
//	return res;
//}

//double gslVelocityDispersionBulge(double z, void* p) {
//	struct gslVelocityDispersionBulgeParams* fp = (struct gslVelocityDispersionBulgeParams*)p;
//	double R2 = gsl_pow_2(fp->R);
//	double z2 = gsl_pow_2(z);
//	double r = gsl_hypot(fp->R, z);
//	double r2 = gsl_pow_2(r);
//	double r3 = gsl_pow_3(r);
//	double haloTemp = fp->densityHalo * gsl_pow_3(fp->rHalo);
//	double hypotbz = gsl_hypot(fp->bDisk, z);
//	double diskTemp = fp->aDisk + hypotbz;
//	//only bulge density:
//	//double temp = fp->mMassBlackHole / r3 + fp->mMassBulge / (r * gsl_pow_2(fp->aBulge + r));
//	//	 + fp->mMassDisk* diskTemp / (hypotbz*pow(R2 + gsl_pow_2(diskTemp), 1.5))
//	//	 -4 * M_PI * haloTemp / (r2 * (fp->rHalo + r)) +4*M_PI*haloTemp *log((fp->rHalo + r) / fp->rHalo)/ r3;
//	//temp = temp* fp->G* z / (r* gsl_pow_3(fp->aBulge + r));
//	//all densities:
//	double temp = fp->aBulge * fp->mMassBulge / (2 * M_PI * r * gsl_pow_3(fp->aBulge + r)) + haloTemp / (r * gsl_pow_2(fp->rHalo + r)) + gsl_pow_2(fp->bDisk) * fp->mMassDisk * (fp->aDisk * R2 + gsl_pow_2(diskTemp) * (fp->aDisk + 3 * hypotbz)) / (4 * M_PI * gsl_pow_2(hypotbz) * pow(R2 + gsl_pow_2(diskTemp), 2.5));
//	temp = temp * fp->G * z * (fp->mMassBlackHole / r + fp->mMassBulge / (r * gsl_pow_2(fp->aBulge + r)) - 4 * M_PI * haloTemp / (r2 * (fp->rHalo + r)) + fp->mMassDisk*diskTemp / (hypotbz*pow(R2 + gsl_pow_2(diskTemp), 1.5)) + 4 * M_PI * haloTemp * log((fp->rHalo + r) / fp->rHalo) / r3);
//	return temp;
//}

//Vec3D MWPotential::sampleDisk(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax){
//	double smallestx = closestToZero(xMin, xMax);
//	double smallesty = closestToZero(yMin, yMax);
//	double smallestz = closestToZero(zMin, zMax);
//
//	double k = densityDisk(smallestx, smallesty, smallestz);
//	std::random_device rd;  //Will be used to obtain a seed for the random number engine
//	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//	std::uniform_real_distribution<> disx(xMin, xMax);
//	std::uniform_real_distribution<> disy(yMin, yMax);
//	std::uniform_real_distribution<> disz(zMin, zMax);
//	std::uniform_real_distribution<> disaccept(0, k);
//
//	while (true) {
//		double x = disx(gen);
//		double y = disy(gen);
//		double z = disz(gen);
//
//		double accept = disaccept(gen);
//		double temp = densityDisk(x, y, z);
//
//		if (accept < temp)
//			return Vec3D(x, y, z);
//	}
//	return Vec3D();
//}

//double MWPotential::gslDensityDiskxNew(double x, void* p) {
//	gsl_function F;
//	F.function = &gslDensityDiskyNew;
//	struct gslDensityConeParamsNew* fp = (struct gslDensityConeParamsNew*)p;
//	fp->x = x;
//	F.params = fp;
//	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
//	double result, error;
//	double temp = fp->z * fp->maxR / fp->maxZ;
//	double boundary = sqrt(temp * temp - x * x);
//	gsl_integration_qag(&F, -boundary, boundary, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
//	gsl_integration_workspace_free(iw);
//	return result;
//}
//
//double MWPotential::gslDensityDiskyNew(double y, void* p) {
//	struct gslDensityConeParamsNew* fp = (struct gslDensityConeParamsNew*)p;
//	Vec3D location = Vec3D(fp->x, y, fp->z);
//	location = *(fp->transformation) * location;
//
//	double z2b2 = location.z * location.z + bDisk * bDisk;
//	double sz2b2 = sqrt(z2b2);
//	double R = gsl_hypot(location.x, location.y);
//	double temp = gsl_pow_2(bDisk) * mMassDisk / (4 * M_PI) * (aDisk * pow(R, 2) + (aDisk + 3 * sz2b2) * pow(aDisk + sz2b2, 2)) / (pow(pow(R, 2) + pow(aDisk + sz2b2, 2), 2.5) * pow(z2b2, 1.5));
//
//	return temp;
//}
//
//double MWPotential::gslDensityDiskzNew(double z, void* p) {
//	gsl_function F;
//	F.function = &gslDensityDiskxNew;
//	struct gslDensityConeParamsNew* fp = (struct gslDensityConeParamsNew*)p;
//	fp->z = z;
//	F.params = fp;
//	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
//	double result, error;
//	double boundary = fp->z * fp->maxR / fp->maxZ;
//	gsl_integration_qag(&F, -boundary, boundary, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
//	gsl_integration_workspace_free(iw);
//	return result;
//
//}

//double MWPotential::massDisk(Matrix* transformation, double minDist, double maxDist, double maxR) {
//	gsl_function F;
//	F.function = &gslDensityDiskzNew;
//
//	gslDensityConeParamsNew densityDiskConeParam{ transformation, maxR, maxDist, 0,0 };
//	F.params = &densityDiskConeParam;
//	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
//	double result, error;
//
//	gsl_integration_qag(&F, minDist, maxDist, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
//	gsl_integration_workspace_free(iw);
//	return result;
//}

//struct gslDensityConeParamsNew { Matrix* transformation; double maxR; double maxZ; double z; double x; };