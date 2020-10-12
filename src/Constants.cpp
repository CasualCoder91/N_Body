#include "Constants.h"

std::string Constants::filePath = "./simulation.cfg";
Configuration config = Configuration(Constants::filePath);

double Constants::softening = config.GetDouble("softening");
double Constants::precission = config.GetDouble("precission");
std::string Constants::title = config.GetString("title");
double Constants::boxLength = config.GetDouble("boxLength");
double Constants::dt = config.GetDouble("dt"); //[day]
int Constants::nTimesteps = config.GetInt("nTimesteps");
int Constants::outputTimestep = config.GetInt("outputTimestep");
int Constants::simulationID = 1;

//cluster mass
double Constants::minMass = config.GetDouble("minMass");
double Constants::maxMass = config.GetDouble("maxMass");
double Constants::alpha = config.GetDouble("alpha");
std::vector<double> Constants::massLimits; //broken powerlaw
std::vector<double> Constants::exponents; //broken powerlaw

//View cone
double Constants::angleOfView = config.GetDouble("angle"); //rad
double Constants::dx = config.GetDouble("dx"); //pc
double Constants::distance = config.GetDouble("distance"); //pc
Vec3D Constants::viewPoint = config.GetVec3D("viewPoint");
Vec3D Constants::focus = config.GetVec3D("focus");

//General Parameters, todo: put into separate class
/** @brief Gravitational constant in astronomical units: parsec/solar mass*km^2/s^2*/
double Constants::G = config.GetDouble("G");// = 4.483e-3;
int Constants::nStars = config.GetInt("nStars");
Vec3D Constants::clusterLocation = config.GetVec3D("offset");

//Analysis Parameters
bool Constants::bEnergy = config.GetBool("bEnergy");
bool Constants::bAverageVelocity = config.GetBool("bAverageVelocity");
bool Constants::bAverage2DVelocity = config.GetBool("bAverage2DVelocity");