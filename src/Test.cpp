#include "..\include\Test.h"

std::string Test::path = "src/Test/";
std::string Test::absolutePath = std::filesystem::current_path().string() +"/" + path;

void Test::pythonScript(std::string fileName){
	std::string filename = path + fileName;
	std::string command = "start python ";
	command += filename;
	system(command.c_str());
}

Test::Test(){}

void Test::sampleFieldStarPositions(int nStars){
	std::cout << "Testing field star positions .." << std::endl;
	sampleFieldStarPositionsOutput(absolutePath,nStars);
	Plot plot = Plot(path, path, true);
	plot.plot("potentialPositions", { std::to_string(nStars) });
}

void Test::sampleFieldStarPositionsOutput(std::string path, int nStars) {
	int firstID = 0;
	std::vector<Star> diskStars = InitialConditions::initStars(firstID, nStars);
	std::vector<Star> bulgeStars = InitialConditions::initStars(firstID, nStars);

	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	MWPotential potential = MWPotential();
	InitialConditions initialConditions = InitialConditions(&potential);
	initialConditions.sampleDiskPositions(diskStars, Vec3D(-40000, -40000, -20000), Vec3D(80000, 80000, 40000));
	initialConditions.sampleBulgePositions(bulgeStars, Vec3D(-6000, -6000, -6000), Vec3D(12000, 12000, 12000));

	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	//positionsDisk.push_back(MWPotential::sampleDisk(-40, 40, -40, 40, -20, 20));
	//positionsBulge.push_back(MWPotential::sampleBuldge(-6, 6, -6, 6, -6, 6));
	InOut::write(diskStars, path + "potentialDiskPositionsSample"+std::to_string(nStars)+".dat");
	InOut::write(bulgeStars, path + "potentialBulgePositionsSample" + std::to_string(nStars) + ".dat");
}

void Test::potentialCircularVelocity() {
	std::cout << "Testing circular velocity .." << std::endl;
	potentialCircularVelocityOutput(path);
	std::cout << " .. plotting data" << std::endl;
	Plot plot = Plot(path, path, true);
	plot.plot("potentialCircularVelocity", {});
	std::cout << " .. done" << std::endl;
}

void Test::potentialCircularVelocityOutput(std::string path){
	std::vector<double> positions;
	std::vector<double> velocities;
	std::vector<double> disk;
	std::vector<double> blackHole;
	std::vector<double> buldge;
	std::vector<double> halo;

	for (double i = 10; i < 30000; i += 100) {
		positions.push_back(i);
		Vec3D position = Vec3D(i, 0, 0);
		velocities.push_back(potential.circularVelocity(&position));
		disk.push_back(potential.circularVelocityDisk(&position));
		blackHole.push_back(potential.circularVelocityBlackHole(&position));
		buldge.push_back(potential.bulgePotential.circularVelocity(&position));
		halo.push_back(potential.circularVelocityHalo(&position));
	}
	std::cout << " .. generating Output" << std::endl;
	InOut::write(positions, velocities, path + "potentialVelocity.dat");
	InOut::write(positions, disk, path + "potentialVelocityDisk.dat");
	InOut::write(positions, blackHole, path + "potentialVelocityBlackHole.dat");
	InOut::write(positions, buldge, path + "potentialVelocityBuldge.dat");
	InOut::write(positions, halo, path + "potentialVelocityHalo.dat");
}

//void Test::initialConditionsMassSalpeterOutput(int nStars) {
//	Constants::nStars = nStars;
//	MWPotential potential = MWPotential();
//	InitialConditions initialConditions = InitialConditions(&potential);
//	int starID = 0;
//	std::vector<Star> stars = initialConditions.initStars(starID, Constants::nStars);
//	initialConditions.initialMassSalpeter(stars, 0.08, 100);
//	std::vector<double> index;
//	std::vector<double> mass;
//	double i = 0;
//	for (const Star& star : stars) {
//		index.push_back(i);
//		i++;
//		mass.push_back(star.mass);
//	}
//
//	InOut::write(index, mass, "initialConditionsMassSalpeter"+std::to_string(nStars)+".dat");
//}
//
//void Test::initialConditionsMassBulgeOutput(double totalMass){
//	MWPotential potential = MWPotential();
//	InitialConditions initialConditions = InitialConditions(&potential);
//	int starID = 0;
//	std::vector<Star> stars = initialConditions.bulgeIMF(totalMass, starID);
//	std::vector<double> index;
//	std::vector<double> mass;
//	double i = 0;
//	for (const Star& star : stars) {
//		index.push_back(i);
//		i++;
//		mass.push_back(star.mass);
//	}
//	InOut::write(index, mass, "initialConditionsMassBulge" + std::to_string(totalMass) + ".dat");
//}
//
//void Test::potentialSurfaceDensityBulge(){
//	MWPotential potential = MWPotential();
//	std::vector<double> surfaceDensity;
//	std::vector<double> radius;
//	for (double R = 100; R < 30000; R = R + 100) {
//		radius.push_back(R);
//		surfaceDensity.push_back(potential.bulgePotential.surfaceDensity(R));
//	}
//	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityBulge.dat");
//}
//
//void Test::potentialSurfaceDensityDisk(){
//	MWPotential potential = MWPotential();
//	std::vector<double> surfaceDensity;
//	std::vector<double> radius;
//	for (double R = 100; R < 30000; R = R + 100) {
//		radius.push_back(R);
//		surfaceDensity.push_back(potential.surfaceDensityDisk(R));
//	}
//	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityDisk.dat");
//}
//
//void Test::initialConditionsSampleDisk(){
//	MWPotential potential = MWPotential();
//	InitialConditions initialConditions = InitialConditions(&potential);
//	double gridResolution = 0.001;
//	Vec3D position = Vec3D(5, 5, 0);
//	Vec3D volumeElement = Vec3D(gridResolution, gridResolution, gridResolution);
//	double massInCell = potential.massDisk(position, volumeElement);
//	int starID = 0;
//	std::vector<Star> starsInCell = initialConditions.diskIMF(massInCell,starID); //stars with mass
//	initialConditions.sampleDiskPositions(starsInCell, position, volumeElement);
//	initialConditions.sampleDiskVelocities(starsInCell);
//}
//
void Test::massDistribution(double z, double dx){
	std::cout << "Testing mass distribution: \n- generating Output " << std::endl;
	massDistributionDiskOutput(path, z, dx); // [z] = [kpc]
	massDistributionBulgeOutput(path, z, dx); // [z] = [kpc]
	std::cout << "- plotting data" << std::endl;
	Plot plot = Plot(path, path, true);
	plot.plot("massDistribution", { std::to_string((int)z) });
	std::cout << "- done" << std::endl;
}

void Test::massDistributionDiskOutput(std::string path, double z, double dx){
	MWPotential potential = MWPotential();
	std::vector<Vec3D> Output;
	double border = 30000;
	ProgressBar progressBar = ProgressBar(-border, border -dx);
	std::cout << "Generating disk mass distribution" << std::endl;
	for (double x = -border; x < border; x += dx) {
		for (double y = -border; y < border; y += dx) {
			double starMass = potential.massDisk(Vec3D(x, y, z), Vec3D(dx, dx, 100));
			Output.push_back(Vec3D(x, y, starMass));
		}
		progressBar.Update(x);
		progressBar.Print();
	}
	std::cout << "done" << std::endl;

	InOut::write(Output, path+"massDistributionDisk_z" + std::to_string((int)z) + ".dat");
}

void Test::massDistributionBulgeOutput(std::string path, double z, double dx){
	MWPotential potential = MWPotential();
	std::vector<Vec3D> Output;
	double border = 30000;
	ProgressBar progressBar = ProgressBar(-border, border-dx);
	std::cout << "Generating disk mass distribution" << std::endl;
	for (double x = -border; x < border; x += dx) {
		for (double y = -border; y < border; y += dx) {
			double starMass = potential.bulgePotential.mass(Vec3D(x, y, z), Vec3D(dx, dx, 100));
			Output.push_back(Vec3D(x, y, starMass));
		}
		progressBar.Update(x);
		progressBar.Print();
	}

	InOut::write(Output, path+"massDistributionBulge_z" + std::to_string((int)z) + ".dat");

}
//
//void Test::massDistributionTimer(){
//	double totalMass = 0;
//	MWPotential potential = MWPotential();
//	for (double x = 0; x < 10000; x = x + 100) {
//		Vec3D corner = Vec3D(x, 0, -50);
//		Vec3D volumeElement = Vec3D(10, 10, 10);
//		Vec3D center = Vec3D::center(corner, corner+volumeElement);
//
//		std::cout << "Volume Element at (center): " << center.length() * 1e-3 << "kpc" << std::endl;
//		std::cout << "Volume Element dx:" << volumeElement.x * 1e-3 << "kpc" << std::endl;
//		//std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
//		double starMass = potential.massDisk(corner, volumeElement)+potential.bulgePotential.mass(corner, volumeElement);
//		totalMass += starMass;
//
//		std::cout << "starMass: " << starMass << std::endl;
//
//		//std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
//		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
//		//std::cout << "Mass inside Volume: " << starMass << std::endl;
//		//std::cout << "Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
//		//startTime = std::chrono::steady_clock::now();
//		//std::vector<Star*> stars = InitialConditions::massDisk(1e4);
//		//endTime = std::chrono::steady_clock::now();
//		//time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
//		//std::cout << "Sample mass Disk - Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
//	}
//	std::cout << "totalMass: " << totalMass << std::endl;
//}
//
//void Test::velocityDispersionBulgerGC(){
//	MWPotential potential = MWPotential();
//	for (double r = 100; r < 6000; r += 100) {
//		std::cout << "r: " << r << " | velocityDistribution: " << potential.velocityDispersionBulge(r) << std::endl;
//	}
//}
//
//void Test::velocityBulge(){
//	InitialConditions conditions = InitialConditions(&potential);
//	std::cout << "Test: velocityBulge() start" << std::endl;
//	double rBulge = 2e3; //pc
//	double dEarth = 8.3e3; //pc
//	int starID = 0;
//	std::vector<double> longitude;
//	std::vector<double> velocityDispersion;
//	std::vector<double> averageVelocity;
//	std::vector<double> bValues{ -4,-6,-8 };
//	for (double b : bValues) {
//		for (double l = -10; l <= 10; l += 0.1) {
//			Vec3D pHGP = Vec3D(1, l * Constants::degInRad, b * Constants::degInRad); //position in Heliocentric Galactic Polar
//			Vec3D vHGP, pHCA, vHCA, pLSR, vLSR,pGCA,vGCA;
//
//			Projection::HGPtoHCA(pHGP, vHGP, pHCA, vHCA);
//			Projection::HCAtoLSR(pHCA, vHCA, pLSR, vLSR);
//			Projection::LSRtoGCA(pLSR, vLSR, pGCA, vGCA);
//			Vec3D focus = pGCA;
//
//			//double x = dEarth - dEarth * cos(b * Constants::degInRad) * cos(l * Constants::degInRad);
//			//double y = dEarth * cos(b * Constants::degInRad) * sin(l * Constants::degInRad);
//			//double z = dEarth * sin(b * Constants::degInRad);
//			//Vec3D focus = Vec3D(x, y, z);
//			double r = focus.length();
//			std::cout << "focus:" << focus.print() << " r:" << r << std::endl;
//			std::vector<Star*> stars = this->initBulgeStars(starID,focus, Vec3D(8300, 0, 27), 7.5e3,0.05);
//			std::cout << "NStars: " << stars.size() << std::endl;
//			std::vector<double> radialVelocities;
//			for (Star* star : stars) {
//				Projection::GCAtoLSR(star->position, star->velocity, pLSR, vLSR);
//				Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);
//				Projection::HCAtoHGP(pHCA, vHCA, pHGP,vHGP);
//				//star->position.x += dEarth; //transform to heleocentric
//				//Vec3D sphericalVelocity = star->velocity.cartesianToSphericalV(star->position);
//				radialVelocities.push_back(vHGP.x);
//			}
//			velocityDispersion.push_back(Analysis::dispersion(radialVelocities));
//			averageVelocity.push_back(Analysis::average(radialVelocities));
//			std::cout << "b:" << b << " l:" << l << " disp:" << velocityDispersion.back()<<" mean:" << averageVelocity.back() << std::endl;
//			for (Star* star : stars) {
//				delete star;
//			}
//			longitude.push_back(l);
//		}
//
//		InOut::write(longitude, velocityDispersion, path + "bulgeDispersion" + std::to_string((int)b) + ".dat");
//		InOut::write(longitude, averageVelocity, path + "bulgeMeanVelocity" + std::to_string((int)b) + ".dat");
//		longitude.clear();
//		velocityDispersion.clear();
//		averageVelocity.clear();
//	}
//	Plot plot = Plot(absolutePath, absolutePath, true);
//	plot.plot("bulgeDispersion", { });
//	plot.plot("bulgeMeanVelocity", { });
//}
//
//void Test::velocityDisk(){
//	static const double Q = 2.7; //Toomre
//	static const double R = 8300; // location of sun
//	static const double surfaceDensity = potential.surfaceDensityDisk(R);
//	static const double freq = potential.epicyclicFrequency(R, 27);
//	static const double k = Q * 3.36 * Constants::G * surfaceDensity / freq * exp(R / potential.aDisk);
//	std::cout << "surfaceDensity: " << surfaceDensity << std::endl;
//	std::cout << "frequency: " << freq << std::endl;
//	std::cout << "constant: " << k << std::endl;
//	printf_s("Dispersion in solar neighborhood: %f\n", k * exp(-R / potential.aDisk));
//	Vec3D position = Constants::positionSun;
//	Vec3D velocity;
//	initialConditions.sampleDiskVelocity(velocity, position);
//	std::cout << "Sampled velocity: " << velocity.print() << std::endl << std::endl;
//}
//
//void Test::bulgeMass(){
//	MWPotential potential = MWPotential();
//	InitialConditions conditions = InitialConditions(&potential);
//	//std::cout << "Disk mass inside aBulge: " << potential.massDisk(Vec3D(-potential.aBulge, -potential.aBulge, -potential.aBulge), Vec3D(2 * potential.aBulge, 2 * potential.aBulge, 2 * potential.aBulge));
//	std::cout << "DensityProfileBulge: start" << std::endl;
//	double rBulge = 4e3; //pc
//	double dEarth = 8.5e3; //pc
//	int starID = 0;
//	std::vector<double> longitude;
//	std::vector<double> mass;
//	std::vector<double> bValues{ 0, 1, 2,3,4,5,6,7,8,9  };
//	//const double degInRad = 0.0174533;
//	for (double b : bValues) {
//		for (double l = -10; l <= 10; l += 1) {
//			longitude.push_back(l);
//			double x = dEarth - dEarth *cos(b * Constants::degInRad) * cos(l * Constants::degInRad);
//			double y = dEarth*cos(b * Constants::degInRad) * sin(l * Constants::degInRad);
//			double z = dEarth*sin(b * Constants::degInRad);
//			Vec3D focus = Vec3D(x, y, z);
//			double r = focus.length();
//			double totalMass = conditions.bulgeStarMass(focus, Vec3D(dEarth,0,0), dEarth+rBulge, 100, 1);
//			mass.push_back(totalMass);
//		}
//		InOut::write(longitude, mass, "testBulgeMass" + std::to_string((int)b) + ".dat");
//		longitude.clear();
//		mass.clear();
//	}
//	Plot plot = Plot(path, path, true);
//	plot.plot("bulgeDispersion", { });
//}
//
//void Test::checkBrokenPowerLaw(){
//	std::cout << "Test: brokenPowerLaw() start" << std::endl;
//	int starID = 0, nStars = 10000;
//	std::vector<Star> stars = this->initialConditions.initStars(starID, nStars);
//	std::vector<double> massLimits = { 0.01,0.08,0.5,1,125 };
//	std::vector<double> exponents = { 0.3,2.0,0.3,2.3 };
//	initialConditions.brokenPowerLaw(stars, massLimits, exponents);
//	std::vector<double> starMass;
//	std::vector<double> index;
//	for (int i = 0; i < stars.size();i++) {
//		starMass.push_back(stars[i].mass);
//		index.push_back(i);
//	}
//	InOut::write(index, starMass, path + "brokenPowerLaw"+std::to_string(nStars)+".dat");
//	std::cout << " .. plotting data" << std::endl;
//	Plot plot = Plot(path, path, true);
//	plot.plot("initialConditionsBrokenPowerLaw", { std::to_string(nStars)});
//	std::cout << " .. done" << std::endl;
//}
//
//void Test::transformation(){
//	MWPotential potential = MWPotential();
//	InitialConditions conditions = InitialConditions(&potential);
//
//	Vec3D rotationVec = (Constants::focus - Constants::viewPoint).normalize();;
//	Vec3D translationVec = Constants::viewPoint;
//
//	Matrix rotationM = Matrix::transformation(rotationVec, translationVec);
//
//	double angleOfView = Constants::angleOfView;
//	double r = 0.5 * Constants::distance * tan(angleOfView);
//	Vec3D boundary1 = Vec3D(-r, -r, 0);
//	Vec3D boundary2 = Vec3D(r,	r, Constants::distance);
//
//	boundary1 = rotationM * boundary1;
//	boundary2 = rotationM * boundary2;
//
//	std::cout << "boundary1: " << boundary1.print() << std::endl;
//	std::cout << "boundary2: " << boundary2.print() << std::endl;
//
//	std::cout << "DiskMass: " << potential.massDisk(&rotationM, 1000, 500) << std::endl;
//	std::cout << "DiskMass: " << potential.massDisk(&rotationM, 1, 2) << std::endl;
//	std::cout << "DiskMass: " << potential.massDisk(&rotationM, 2, 1) << std::endl;
//	std::cout << "DiskMass: " << potential.massDisk(&rotationM, 3, 1) << std::endl;
//	std::cout << "DiskMass: " << potential.massDisk(&rotationM, 4, 1) << std::endl;
//}
//
////void Test::velocityBulge(){
////	Parameters parameters = Parameters();
////	MWPotential potential = MWPotential(&parameters);
////	std::cout << "BulgeVelocityTest: start" << std::endl;
////	std::vector<double> longitude;
////	std::vector<double> dispersion;
////	std::vector<double> bValues{ -4, -6, -8 };
////	for (double b : bValues) {
////		for (double l = -10; l < 10; l += 0.2) {
////			longitude.push_back(l);
////			std::vector<Vec3D*> velocities;
////			double distance = 8490;
////			double x = (8500 - distance * cos(b * 0.0174533) * cos(l * 0.0174533));
////			double y = distance * cos(b * 0.0174533) * sin(l * 0.0174533);
////			double z = distance * sin(b * 0.0174533);
////			double r = gsl_hypot3(x, y, z);
////			double R = gsl_hypot(x, y);
////			double disp = potential.velocityDispersionBulge(r);
////			//std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
////			dispersion.push_back(disp);
////			//std::cout << "radial dist:" << MWPotential::radialVelocityDispersionBulge(R, z) << std::endl;
////		}
////		InOut::write(longitude, dispersion, path + "bulgeDispersion" + std::to_string((int)b)+".dat");
////		longitude.clear();
////		dispersion.clear();
////	}
////	Plot plot = Plot(absolutePath, absolutePath, true);
////	plot.plot("bulgeDispersion", { });
////}
//
//void Test::velocityBulgeR() {
//	std::vector<double> radius;
//	std::vector<double> dispersion;
//	MWPotential potential = MWPotential();
//	for (double r = 10; r < 1000; r += 1) {
//		radius.push_back(r);
//		double disp = potential.velocityDispersionBulge(r);
//		std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
//		dispersion.push_back(disp);
//		//std::cout << "radial dist:" << MWPotential::radialVelocityDispersionBulge(R, z) << std::endl;
//	}
//	InOut::write(radius, dispersion, "bulgeDispersion.dat");
//
//}
//
//void Test::initialConditionsSampleBulgeVelocity(){
//	Constants::nStars = 1;
//	MWPotential potential = MWPotential();
//	InitialConditions initialConditions = InitialConditions(&potential);
//	int starID = 0;
//	std::vector<Star> stars = initialConditions.initStars(starID, Constants::nStars);
//	double delta = 100;
//	std::vector<double> bValues{ -4, -6, -8 };
//	std::vector<double> longitude;
//	std::vector<double> averageVelocity;
//	for (double b : bValues) {
//		for (double l = -10; l < 10; l += 0.2) {
//			longitude.push_back(l);
//			double distance = 8490;
//			double x = (8500 - distance * cos(b * 0.0174533) * cos(l * 0.0174533));
//			double y = distance * cos(b * 0.0174533) * sin(l * 0.0174533);
//			double z = distance * sin(b * 0.0174533);
//			double r = gsl_hypot3(x, y, z);
//			initialConditions.sampleBulgePositions(stars, Vec3D(r, 0, 0), Vec3D(delta, delta, delta));
//			initialConditions.sampleBulgeVelocities(stars);
//			std::vector<Vec3D> velocities;
//			for (const Star& star : stars) {
//				velocities.push_back(star.velocity);
//			}
//			//std::cout << "r: " << r << " | average velocity: " << Analysis::average(velocities) << " | Dispersion: " << Analysis::dispersion(velocities) << std::endl;
//			averageVelocity.push_back(Analysis::average(velocities));
//		}
//		InOut::write(longitude, averageVelocity, path + "bulgeMeanVelocity" + std::to_string((int)b) + ".dat");
//		longitude.clear();
//		averageVelocity.clear();
//	}
//	Plot plot = Plot(path, path, true);
//	plot.plot("bulgeMeanVelocity", { });
//}
//
//void Test::escapeVelocity(){
//	double delta = 50;
//	MWPotential potential = MWPotential();
//	for(double x = 0; x < 10000; x += delta) {
//		std::cout << "r:" << x << " | velocity: " << potential.escapeVelocity(&Vec3D(x, 0, 0)) << std::endl;
//	}
//}
//
//void Test::initialConditionsInitFieldStars(){
//	std::cout << "Testing field star creation" << std::endl;
//	MWPotential potential = MWPotential();
//	InitialConditions initialConditions = InitialConditions(&potential);
//	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
//	int starID = 0;
//	std::vector<Star> stars = initialConditions.initFieldStars(starID,Constants::focus, Constants::viewPoint, Constants::distance, Constants::angleOfView);
//	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
//	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
//	std::cout << "Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
//	InOut::write(stars, path+"fieldStars.dat");
//}
//
//std::vector<Star*> Test::initBulgeStars(int& starID, Vec3D focus, Vec3D viewPoint, double distance, double angleOfView){
//	std::cout << "Test::initBulgeStars() start" << std::endl;
//	Vec3D direction = (focus - viewPoint).normalize();
//	double angleOfViewRad = angleOfView * Constants::degInRad; //convert degrees in rad
//	double coneR = distance * tan(0.5 * angleOfViewRad);
//	Matrix transformationMatrix = Matrix::transformation(direction, viewPoint);
//	Vec3D coneBoundaryMin = transformationMatrix * Vec3D(-coneR, -coneR, 0);
//	Vec3D coneBoundaryMax = transformationMatrix * Vec3D(coneR, coneR, distance * 1.01); //1.01 to make sure boundary is not inside cone.
//	initialConditions.setBoundaries(coneBoundaryMin, coneBoundaryMax);
//
//	std::vector<Star*> fieldStars; //return vector
//
//	double bulgeMass = potential.bulgePotential.mass(&transformationMatrix, distance, coneR);
//	std::vector<Star> bulgeStars = initialConditions.bulgeIMF(bulgeMass, starID);
//	if (bulgeStars.size() > 0) {
//		//sampleWang(bulgeStars, corner, volumeElement);
//		initialConditions.sampleBulgePositions(bulgeStars, coneBoundaryMin, coneBoundaryMax, coneR, distance, &transformationMatrix); //test
//		initialConditions.sampleBulgeVelocities(bulgeStars);
//		fieldStars.insert(std::end(fieldStars), std::begin(bulgeStars), std::end(bulgeStars));
//	}
//	std::cout << "Test::initBulgeStars() done" << std::endl;
//	return fieldStars;
//}
