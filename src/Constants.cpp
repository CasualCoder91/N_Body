#include "Constants.h"

std::string Constants::filePath = "./simulation.cfg";
Configuration config = Configuration(Constants::filePath);

double Constants::softening = config.GetDouble("softening");
double Constants::precission = config.GetDouble("precission");
std::string Constants::title = config.GetString("title");
double Constants::plummer_radius = config.GetDouble("plummer_radius");
double Constants::dt = config.GetDouble("dt"); //[day]
int Constants::nTimesteps = config.GetInt("nTimesteps");
int Constants::outputTimestep = config.GetInt("outputTimestep");

//cluster
std::vector<double> Constants::massLimits = config.GetDoubleVector("massLimits"); //broken powerlaw
std::vector<double> Constants::exponents = config.GetDoubleVector("exponents"); //broken powerlaw
Vec3D Constants::clusterLocation = config.GetVec3D("offset");
bool Constants::bMcLuster = config.GetBool("bMcLuster");
std::string Constants::mcluster_filepath = config.GetString("mcluster_filepath");

//View cone
double Constants::angleOfView = config.GetDouble("angle"); //degree
double Constants::distance = config.GetDouble("distance"); //pc
Vec3D Constants::viewPoint = config.GetVec3D("viewPoint");
Vec3D Constants::focus = config.GetVec3D("focus");

//General Parameters, todo: put into separate class
/** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
double Constants::G = config.GetDouble("G");// = 4.483e-3;
double Constants::degInRad = 0.017453292519943295769480039544949184171511371924940265706568;
double Constants::radInArcsec = 206264.806;
double Constants::radmyrInArcsecyr = 0.206264806;
double Constants::kmInpc = 3.086e-13;
int Constants::nStars = config.GetInt("nStars");

//Analysis Parameters
//bool Constants::bEnergy = config.GetBool("bEnergy");
//bool Constants::bAverageVelocity = config.GetBool("bAverageVelocity");
//bool Constants::bAverage2DVelocity = config.GetBool("bAverage2DVelocity");

//Transformation
Vec3D Constants::positionSun = Vec3D(8300, 0, 27).cartesianToCylindrical();
Vec3D Constants::velocitySun = Vec3D(11.1, 12.24, 7.25); //kms
double Constants::circularVelocitySun = -233.1;

double Constants::pi = 3.14159265358979323846;
double Constants::pi2 = 2.*3.14159265358979323846;
double Constants::pi4 = 4.*3.14159265358979323846;

double Constants::eps_magnitude = config.GetDouble("eps_magnitude");