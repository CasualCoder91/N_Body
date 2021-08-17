/**
 * Barens Hut Simulation
 *
 * Stars are initialized (mass, possition- and velocity distribution). Octree is initialized based on the stars.
 * Stars accelerations are updated via the octree. Velocity and possition of the stars are updated via an integrator.
 * Parts of the code can be run in parallel via OpenMP
 *
 * @author Alarich Herzner
 * @version 0.9 05.03.2020
*/

#include <iostream>
#include <random>
#include <vector>
#include <chrono> //for timer
#include <direct.h> //_mkdir

#include "Star.h"
#include "Node.h"
#include "InOut.h"
#include "Integrator.h"
#include "InitialConditions.h"
#include "Analysis.h"
#include "Database.h"
#include "Simulation.h"
#include "MWPotential.h"
#include "WangPotential.h"
#include "Projection.h"

#include "Test.h"
#include "Plot.h"

typedef int mode_t;
//namespace plt = matplotlibcpp;

//global parameters
bool debug = false;

int main() {

	//InitialConditions testInitialConditions = InitialConditions();
	//MWPotential testPotential = MWPotential(Vec3D(0, 0, 0));
	//double boxSize = 0.002;
	//std::vector<Star*> stars = testInitialConditions.initDiskStars(0, Vec3D(0.1, boxSize, 0), Vec3D(0.1+boxSize, 0, 0), 0.1, &testPotential);
	//InOut::write(stars, "testInitialConditions.dat");

	//Test::potentialSurfaceDensityBulge();
	//Test::potentialSurfaceDensityDisk();

	//Test::potentialCircularVelocityOutput();
	//Test::velocityBulge();
	//Test::potentialCircularVelocityOutput();
	//Test::initialConditionsSampleBulgeVelocity();
	//Test::escapeVelocity();
	//Test::initialConditionsInitFieldStars();
	//MWPotential::generateVelocityDistributionBulgeLookupTable(25000);
	//std::vector<std::vector<double>> test = InOut::readDoubleMatrix("velocityDistributionBulgeTable.dat");

	Database db = Database();
	db.open();
	db.setup();

	while (true) {
		int selection;
		int simulationID = -1;
		std::cout << "[1] Load Simulation\n[2] New Simulation\n[3] Run Tests\n[4] Exit" << std::endl;
		std::cin >> selection;
		std::cin.clear();
		if (selection == 1) {
			if(!db.printSimulations())
				continue;
			std::cout << "Input ID to select simulation" << std::endl;
			int simulationID = 0;
			std::cin >> simulationID;
			std::cin.clear();
			db.selectConstants(simulationID); //Load config used for specified simulation
			Simulation simulation = Simulation(simulationID,&db);
			std::cout << "[1] Ouput\n[2] Analysis\n[3] Generate observables (HTP/magnitude)" << std::endl;
			std::cin >> selection;
			std::cin.clear();
			if (selection == 1) {
				std::cout << "[1] 3D\n[2] 2D" << std::endl;
				std::cin >> selection;
				std::cin.clear();
				std::string directory = "Simulation" + std::to_string(simulation.getID());
				directory = InOut::makeDirectory(directory);
				std::cout << "Files will be written to: " << directory << std::endl;
				if (selection == 1) {
					db.outputStars(simulation.getID(), directory + "/stars", false, true, true);
				}
				else if (selection == 2) {
					db.outputStars2D(simulation.getID(), directory + "/stars");
				}
				std::cout << "done" << std::endl;
			}
			else if(selection==2) {
				Analysis analysis = Analysis(simulationID,&db);
				std::cout << "What would you like to analyze?" << std::endl;
				std::cout << "[1] Energy\n[2] Velocity\n[3] velocityHTP\n[4] Cluster\n[5] Map observed stars" << std::endl;
				std::cin >> selection;
				std::cin.clear();
				std::vector<int> timeSteps = db.selectTimesteps(simulationID);
				if (selection == 1) {//Energy
					analysis.energy();
					std::cout << "Energy analysis done" << std::endl;
				}
				else if (selection == 2) { //Velocity3D
					for (int timeStep : timeSteps) {
						//std::vector<Star*> stars = db.selectStars(simulationID, timeStep);
						std::vector<Vec3D> clusterVelocities = db.selectVelocities3D(simulationID, timeStep, false, true);
						std::vector<Vec3D> fsVelocities = db.selectVelocities3D(simulationID, timeStep, true, false);
						double avgVel3DCluster = analysis.average(clusterVelocities);
						double disp3DCluster = analysis.dispersion(clusterVelocities, avgVel3DCluster);
						double avgVel3DFS = analysis.average(fsVelocities);
						double disp3DFS = analysis.dispersion(fsVelocities, avgVel3DFS);
						db.insertAnalysisdtVelocity3D(simulationID, timeStep, avgVel3DCluster,disp3DCluster,avgVel3DFS,disp3DFS);
					}
					std::cout << "3D velocity analysis done" << std::endl;
				}
				else if (selection == 3) { //velocityHTP
					for (int timeStep : timeSteps) {
						std::vector<Vec2D> clusterVelocities = db.selectVelocitiesHTP(simulationID, timeStep, false, true, Constants::minMagnitude);
						std::vector<Vec2D> fsVelocities = db.selectVelocitiesHTP(simulationID, timeStep, true, false, Constants::minMagnitude);
						double avgVelHTPCluster = analysis.average(clusterVelocities);
						double dispHTPCluster = analysis.dispersion(clusterVelocities, avgVelHTPCluster);
						double avgVelHTPFS = analysis.average(fsVelocities);
						double dispHTPFS = analysis.dispersion(fsVelocities, avgVelHTPFS);
						db.insertAnalysisdtVelocity2D(simulationID, timeStep, avgVelHTPCluster, dispHTPCluster, avgVelHTPFS, dispHTPFS);
					}
					std::cout << "HTP velocity analysis done" << std::endl;
				}
				else if (selection == 4) {//Cluster
					std::cout << "[1] Simulated stars\n[2] Observed stars" << std::endl;
					std::cin >> selection;
					std::cin.clear();

					std::cout << "Running cluster analysis ..." << std::endl;
					if (selection == 1) {
						analysis.cluster(db.select_points(simulationID, 0, -1));
					}
					else {
						analysis.cluster(db.select_points(simulationID, 0, -1,true));
					}
					std::cout << "Cluster analysis done" << std::endl;
				}
				else if (selection == 5) {
					analysis.map_observed();
					std::cout << "Mapping done!" << std::endl;
				}
				else {
					std::cout << "Feature not yet implemented" << std::endl;
				}
			}
			else {
				std::cout << "[1] Simulated stars\n[2] Observed stars" << std::endl;
				std::cin >> selection;
				std::cin.clear();
				if (selection == 1) {
					std::cout << "generating HTP positions ..." << std::endl;
					db.generateHTP(simulationID,false); //todo: move this to Analysis
					std::cout << "generating magnitude ..." << std::endl;
					db.generateMagnitude(simulationID);
					//HTP velocity
					std::cout << "generating HTP velocities ..." << std::endl;
					Analysis analysis = Analysis(simulationID, &db);
					analysis.generateHTPVelocity(false,false);
					std::cout << "done\n" << std::endl;
				}
				else if (selection == 2) {
					std::cout << "generating HTP velocities ..." << std::endl;
					Analysis analysis = Analysis(simulationID, &db);
					analysis.generateHTPVelocity(true);
					std::cout << "done\n" << std::endl;
				}
			}
		}
		else if (selection == 2) {
			simulationID = db.insertSimulation();
			db.insertAnalysis(simulationID);
			Simulation simulation = Simulation(simulationID, &db);
			std::cout << "New simulation created. ID = " << simulationID << std::endl;
			std::cout << "Starting simulation" << std::endl;
			simulation.run();
		}
		else if (selection == 3) {
			//Vec3D pGCP = Vec3D(-9000, 206.3059*Constants::degInRad, - 02.0720 * Constants::degInRad);
			//Vec3D vGCP = Vec3D();

			//Vec3D pGCA = Vec3D(0.1,0.2,0.3), vGCA = Vec3D(-1,-2,-3);
			//std::cout << "GCA: " << pGCA.print() << " | " << vGCA.print() << std::endl;
			//Vec3D pLSR, vLSR;
			//Projection::GCAtoLSR(pGCA, vGCA, pLSR, vLSR);
			//Projection::LSRtoGCA(pLSR, vLSR, pGCA, vGCA);
			//std::cout << "GCA: " << pGCA.print() << " | " << vGCA.print() << std::endl;

			//Vec3D pGCA = Vec3D(9594, -640, -52);
			//Vec3D vGCA = Vec3D(58.62, -12.39, -14.55);

			////Projection::GCPtoGCA(pGCP, vGCP, pGCA, vGCA);
			//std::cout << "GCA: " << pGCA.print() << " | " << vGCA.print() << std::endl;

			//Vec3D pLSR, vLSR;
			//Projection::GCAtoLSR(pGCA, vGCA, pLSR, vLSR);
			//std::cout << "LSR: " << pLSR.print() << " | " << vLSR.print() << std::endl;

			//Vec3D pHCA, vHCA;
			//Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);
			//std::cout << "HCA: " << pHCA.print() << " | " << vHCA.print() << std::endl;

			//Vec3D pHEQ, vHEQ;
			//Projection::HCAtoHEQ(pHCA, vHCA, pHEQ, vHEQ);
			//std::cout << "HEQ: " << pHEQ.print() << " | " << vHEQ.print() << std::endl;

			//Vec3D pHGP, vHGP;
			//Projection::HCAtoHGP(pHCA, vHCA, pHGP, vHGP);
			//std::cout << "HGP: " << pHGP.print() << " | " << vHGP.print() << std::endl;

			//Projection::HGPtoHEQ(pHGP, vHGP, pHEQ, vHEQ);
			//std::cout << "HEQ: " << pHEQ.print() << " | " << vHEQ.print() << std::endl;


			Vec3D pHGP = Vec3D(1719, 16.94 * Constants::degInRad, 0.8 * Constants::degInRad);
			Vec3D vHGP = Vec3D();

			Vec3D pHCA, vHCA, pLSR, vLSR, pGCA, vGCA, pHEQ, vHEQ, pGCP, vGCP;
			Projection::HGPtoHCA(pHGP, vHGP, pHCA, vHCA);
			Projection::HCAtoHEQ(pHCA, vHCA, pHEQ, vHEQ);
			//pHEQ.y = pHEQ.y / Constants::degInRad;
			//pHEQ.z = pHEQ.z / Constants::degInRad;
			std::cout << "HEQ: " << pHEQ.print() << " | " << vHEQ.print() << std::endl;
			Projection::HCAtoLSR(pHCA, vHCA, pLSR, vLSR);
			Projection::LSRtoGCA(pLSR, vLSR, pGCA, vGCA);
			std::cout << "GCA: " << pGCA.print() << " | " << vGCA.print() << std::endl;
			Projection::GCAtoGCP(pGCA, vGCA, pGCP, vGCP);
			std::cout << "GCP: " << pGCP.print() << " | " << vGCP.print() << std::endl;


			//Test test = Test();
			//test.velocityDisk();
			//test.checkBrokenPowerLaw();
			//test.potentialCircularVelocity();
			//test.transformation();
			//Test::massDistribution(500,15000);
			//Test::sampleFieldStarPositions(200);
			//Test::velocityBulgeR();
			//Test::bulgeMass();
			//test.velocityBulge(); // very time intensive

			////double r11 = cos(a) * cos(d);
			////double r12 = -sin(a);
			////double r13 = cos(a) * sin(d);
			////double r21 = sin(a) * cos(d);
			////double r22 = cos(a);
			////double r23 = sin(a) * sin(d);
			////std::cout << r11 << " " << r12 << " " << r13 << std::endl;
			////std::cout << r21 << " " << r22 << " " << r23 << std::endl;

			//double l = 3.366033268750004;
			//double b = 0.47347728280415174;

			//double x = cos(b) * cos(l);
			//double y = cos(b) * sin(l);
			//double z = sin(b);

			//Vec3D positionNCP = Vec3D(x, y, z).normalize();
			//std::cout << "positionNCP: " << positionNCP.print() << std::endl;

			////std::cout << "mNCP: " << std::endl << mNCP << std::endl;

			//l = 4.6496443937754925;
			//b = -0.5050284723570792;

			//x = cos(b) * cos(l);
			//y = cos(b) * sin(l);
			//z = sin(b);

			//Vec3D positionEqui = Vec3D(x, y, z).normalize();

			//Vec3D cross = Vec3D::crossProduct(&positionNCP, &positionEqui);

			//std::cout << "cross: " << cross.print() << std::endl;

			//std::cout << "positionEqui: " << positionEqui.print() << std::endl;


			//std::cout << "mEqui: " << std::endl << mEqui << std::endl;




			//std::cout << m << std::endl << std::endl;
			//std::cout << m << std::endl;

			////Test::wangPositions();
			////Test::checkBrokenPowerLaw();
			//std::cout << WangPotential::ANLM(1, 0, 0) << std::endl;
			//std::cout << WangPotential::totalMass(-5e3, 5e3) << std::endl;
		}
		else if (selection == 4) {
			return 0;
		}
	}
	return 0;
}