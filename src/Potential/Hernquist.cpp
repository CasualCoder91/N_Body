#include "Potential/Hernquist.h"

double Hernquist::mMass = 1.8e10; // SolarMassUnit
double Hernquist::mScaleLength = 0.9e3; // pc

struct gslRParam { double R; };
struct gslDensityConeParams { Matrix* transformation; double distance; double r; double x; double y; };

Hernquist::Hernquist(double mass, double scaleLength){
	mMass = mass;
	mScaleLength = scaleLength;
}

double Hernquist::density(double r){
	return mMass / (2 * M_PI) * mScaleLength / (r * gsl_pow_3(r + mScaleLength));
}

double Hernquist::density(double R, double z){
	return density(gsl_hypot(R, z));
}

double Hernquist::density(double x, double y, double z){
	return density(gsl_hypot3(x,y,z));
}

double Hernquist::gslDensity(double x[], size_t dim, void* p){
	if (dim != 3) {
		fprintf(stderr, "error: dim != 3");
		abort();
	}
	return density(x[0], x[1], x[2]);
}

double Hernquist::gslDensity(double z, void* p){
	struct gslRParam* fp = (struct gslRParam*)p;
	return density(fp->R,z);
}

double Hernquist::potential(double r){
	return - Constants::G*mMass /(r+ mScaleLength);
}

double Hernquist::forceTemp(double r){
	return mMass/(r*gsl_pow_2(mScaleLength+r));
}

double Hernquist::gslDensityX(double x, void* p) {
	gsl_function F;
	F.function = &gslDensityY;
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

double Hernquist::gslDensityY(double y, void* p) {
	gsl_function F;
	F.function = &gslDensityZ;
	struct gslDensityConeParams* fp = (struct gslDensityConeParams*)p;
	fp->y = y;
	F.params = fp;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_integration_qag(&F, fp->distance / fp->r * sqrt(fp->x * fp->x + y * y), fp->distance, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return result;
}

double Hernquist::gslDensityZ(double z, void* p){
	struct gslDensityConeParams* fp = (struct gslDensityConeParams*)p;
	Vec3D location = Vec3D(fp->x, fp->y, z);
	location = *(fp->transformation) * location;
	return density(location.x, location.y, location.z);
}

double Hernquist::mass(Vec3D position, Vec3D volumeElement){
	gsl_monte_function F;

	F.f = &gslDensity;
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
	gsl_rng_free(r);

	return mass;
}

double Hernquist::mass(Matrix* transformationMatrix, double distance, double coneR){
	gsl_function F;
	F.function = &gslDensityX;
	gslDensityConeParams densityDiskConeParam{ transformationMatrix, distance, coneR,0,0 };
	F.params = &densityDiskConeParam;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_integration_qag(&F, -coneR, coneR, 1e-3, 1e-3, 1000, 1, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return result;
}

double Hernquist::surfaceDensity(double R){
	if (R == 0) {
		std::cout << "R=0 not allowed." << std::endl;
		return 0;
	}

	gsl_function F;
	F.function = &gslDensity;
	gslRParam params = { R };
	F.params = &params;
	gsl_integration_workspace* iw = gsl_integration_workspace_alloc(1000);
	double result, error;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, iw, &result, &error);
	gsl_integration_workspace_free(iw);
	return 2.0 * result;
}

double Hernquist::circularVelocity(const Vec3D* position){
	double r = position->length();
	return sqrt(Constants::G * mMass * r) / (r + mScaleLength);
}

double Hernquist::escapeVelocity(double r){
	return sqrt(-2. * potential(r));
}

double Hernquist::escapeVelocity(Vec3D* position){
	return escapeVelocity(position->length());
}

double Hernquist::potentialdr(double r){
	return mMass / gsl_pow_2(mScaleLength + r);
}

double Hernquist::potentialdr2(double r){
	return 2.*mMass / gsl_pow_3(mScaleLength + r);
}
