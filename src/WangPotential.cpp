#include "..\include\WangPotential.h"

const double WangPotential::density0 = 2.1242;
const double WangPotential::q = 0.6;
const double WangPotential::z0 = 0.4e3; //pc
const double WangPotential::x0 = 1.49e3; //pc
const double WangPotential::y0 = 0.58e3; //pc
const double WangPotential::rs = 1e3; //pc
const double WangPotential::mass = 2e10;

struct gslExpansionCoeffParams { unsigned int n; unsigned int l; unsigned int m; };


double WangPotential::density(double xParam, double yParam, double zParam){
	//the density is rotated around z axis by 20 degrees
	double angle = 0.349066; //20 degree
	double x = xParam* cos(angle) - yParam * sin(angle);
	double y = xParam* sin(angle) + yParam * cos(angle);
	double z = zParam;
	double r1 = pow(gsl_pow_2(gsl_pow_2(x / x0) + gsl_pow_2(y / y0)) + gsl_pow_4(z / z0), 0.25);
	double r2 = sqrt((gsl_pow_2(q) * (gsl_pow_2(x) + gsl_pow_2(y)) + gsl_pow_2(z)) / gsl_pow_2(z0));
	if (r1 == 0)//x=y=z=0
		return density0;
	return density0 * (exp(-gsl_pow_2(r1) / 2) + pow(r2, -1.85) * exp(-r2));
}

double WangPotential::gslDensity(double x[], size_t dim, void* p){
	return density(x[0],x[1],x[2]);
}

double WangPotential::expansionCoeffIntegrand(double x[], size_t dim, void* p){
	struct gslExpansionCoeffParams* params = (struct gslExpansionCoeffParams*)p;

	double r = (1 + x[0]) / (1 - x[0]);
	double theta = acos(x[1]);
	double a = r * sin(theta) * cos(x[2]); //x
	double b = r * sin(theta) * sin(x[2]); //y
	double c = r * cos(theta); //z
	double s = r / rs;
	double temp = 2. / gsl_pow_2(1. - x[0]) * WangPotential::density(a, b, c);
	temp = temp * pow(s, params->l) / pow(1 + s, 2 * params->l + 1) * gsl_sf_gegenpoly_n(params->n, 2 * (double)params->l + 1.5, x[0]) * gsl_pow_2(r);
	temp = temp * gsl_sf_legendre_Plm(params->l, params->m, x[1]) * cos(params->m * x[2]);
	return 2. / gsl_pow_2(1. - x[0]) * WangPotential::density(a, b, c)*
		pow(s, params->l) / pow(1 + s, 2 * params->l + 1) * gsl_sf_gegenpoly_n(params->n, 2 * (double)params->l + 1.5, x[0])*gsl_pow_2(r)*
		gsl_sf_legendre_Plm(params->l, params->m, x[1])*cos(params->m*x[2]);
}

double WangPotential::totalMass(double min, double max){
	gsl_monte_function F;
	F.f = &gslDensity;
	F.dim = 3;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double xl[3] = { min, min, min };
	double xu[3] = { max, max, max };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double result, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,
		&result, &err);
	gsl_monte_plain_free(s);
	return result;
}

double WangPotential::PotentialNLM(unsigned int n, unsigned int l, unsigned int m, Vec3D position){
	double s = position.length() / rs;
	return pow(s, l) / pow(1 + s, 2 * l + 1) * gsl_sf_gegenpoly_n(n, 2 * (double)l + 1.5, (s - 1) / (s + 1)) * gsl_sf_legendre_Plm(l, m, cos(position.theta())) * cos(m * position.phi());
}

double WangPotential::distributionFunction(double mass, Vec3D position, Vec3D velocity){
	//double angle = -0.349066; //20 degree
	//double x = position.x * cos(angle) - position.y * sin(angle);
	//double y = position.x * sin(angle) + position.y * cos(angle);
	//position.x = x;
	//position.y = y;
	//double xtemp = velocity.x;
	//double ytemp = velocity.y;
	//velocity.x = xtemp * cos(angle) - ytemp * sin(angle);
	//velocity.y = xtemp * sin(angle) + ytemp * cos(angle);


	Vec3D sigma1 = Vec3D(60, 60, 50); //prograde
	Vec3D sigma2 = Vec3D(60, 60, 50); //retrograde
	Vec3D sigma3 = Vec3D(100, 100, 60); //hot

	Vec3D mC = velocity.cartesianToCylindricalV(position); //momentum
	Vec3D pC = position.cartesianToCylindrical(); //position
	double vcl = 250 / (1 + pow(100 / pC.x, 0.2));

	double temp = exp(-gsl_pow_2(mC.x) / (2 * gsl_pow_2(sigma1.x)) - gsl_pow_2(mC.y - vcl) / (2 * gsl_pow_2(sigma1.y)) - gsl_pow_2(mC.z) / (2 * gsl_pow_2(sigma1.z))) / (pow(2 * M_PI, 1.5) * sigma1.x * sigma1.y * sigma1.z);//prograde
	temp += exp(-gsl_pow_2(mC.x) / (2 * gsl_pow_2(sigma2.x)) - gsl_pow_2(mC.y + vcl) / (2 * gsl_pow_2(sigma2.y)) - gsl_pow_2(mC.z) / (2 * gsl_pow_2(sigma2.z))) / (pow(2 * M_PI, 1.5) * sigma2.x * sigma2.y * sigma2.z);//retrograde
	temp += exp(-gsl_pow_2(velocity.x) / (2 * gsl_pow_2(sigma3.x)) - gsl_pow_2( velocity.y) / (2 * gsl_pow_2(sigma3.y)) - gsl_pow_2( velocity.z) / (2 * gsl_pow_2(sigma3.z))) / (pow(2 * M_PI, 1.5) * sigma3.x * sigma3.y * sigma3.z);//hot
	temp = temp * density(position.x, position.y, position.z);
	return temp;
}

double WangPotential::ANLM(unsigned int n, unsigned int l, unsigned int m){
	double Knl = 0.5 * n * (n + 4 * l + 3) + (l + 1) * (2 * l + 1);
	double kDelta = 0;
	if (m == 0)
		kDelta = 1;
	double Inl = Knl/pow(2.,8.*l+6.)*gsl_sf_gamma(n + 4 * l + 3)/
		(gsl_sf_fact(n) * (n + 2 * l + 1.5) * pow(gsl_sf_gamma(2 * l + 1.5), 2))
		*(1.+kDelta)*M_PI*2./(2.*l+1)* gsl_sf_fact(l+m)/ gsl_sf_fact(l - m);

	gsl_monte_function F;
	gslExpansionCoeffParams params = { n,l,m };
	F.f = &expansionCoeffIntegrand;
	F.dim = 3;
	F.params = &params;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	double xl[3] = { -1, -1, 0 };
	double xu[3] = { 1, 1, 2*M_PI };
	size_t calls = 500000;

	gsl_monte_plain_state* s = gsl_monte_plain_alloc(3);
	double result, err;
	gsl_monte_plain_integrate(&F, xl, xu, 3, calls, r, s,
		&result, &err);
	gsl_monte_plain_free(s);
	return 1 / Inl * result;
}
