#include "..\include\Test.h"

std::string Test::path = "/src/Test/";
std::string Test::absolutePath = std::filesystem::current_path().string() + path;

void Test::pythonScript(std::string fileName){
	std::string filename = path + fileName;
	std::string command = "start python ";
	command += filename;
	system(command.c_str());
}

void Test::sampleFieldStarPositions(int nStars){
	std::cout << "Testing field star positions .." << std::endl;
	sampleFieldStarPositionsOutput(absolutePath,nStars);
	Plot plot = Plot(path, path, true);
	plot.plot("potentialPositions", { std::to_string(nStars) });
}

void Test::sampleFieldStarPositionsOutput(std::string path, int nStars) {
	int firstID = 0;
	std::vector<Star*> diskStars = InitialConditions::initStars(firstID, nStars);
	std::vector<Star*> bulgeStars = InitialConditions::initStars(firstID, nStars);

	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	InitialConditions initialConditions = InitialConditions(&potential);
	initialConditions.sampleDiskPositions(diskStars, Vec3D(-40000, -40000, -20000), Vec3D(80000, 80000, 40000));
	initialConditions.sampleBulgePositions(bulgeStars, Vec3D(-6000, -6000, -6000), Vec3D(12000, 12000, 12000));

	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	//positionsDisk.push_back(Potential::sampleDisk(-40, 40, -40, 40, -20, 20));
	//positionsBulge.push_back(Potential::sampleBuldge(-6, 6, -6, 6, -6, 6));
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
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);

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
		buldge.push_back(potential.circularVelocityBulge(&position));
		halo.push_back(potential.circularVelocityHalo(&position));
	}
	std::cout << " .. generating Output" << std::endl;
	InOut::write(positions, velocities, path + "potentialVelocity.dat");
	InOut::write(positions, disk, path + "potentialVelocityDisk.dat");
	InOut::write(positions, blackHole, path + "potentialVelocityBlackHole.dat");
	InOut::write(positions, buldge, path + "potentialVelocityBuldge.dat");
	InOut::write(positions, halo, path + "potentialVelocityHalo.dat");
}

void Test::initialConditionsMassSalpeterOutput(int nStars) {
	Parameters parameters = Parameters();
	parameters.setNStars(nStars);
	Potential potential = Potential(&parameters);
	InitialConditions initialConditions = InitialConditions(&potential);
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initStars(starID,parameters.getNStars());
	initialConditions.initialMassSalpeter(stars, 0.08, 100);
	std::vector<double> index;
	std::vector<double> mass;
	double i = 0;
	for (Star* star : stars) {
		index.push_back(i);
		i++;
		mass.push_back(star->mass);
	}

	InOut::write(index, mass, "initialConditionsMassSalpeter"+std::to_string(nStars)+".dat");
}

void Test::initialConditionsMassBulgeOutput(double totalMass){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	InitialConditions initialConditions = InitialConditions(&potential);
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initialMassBulge(totalMass, starID);
	std::vector<double> index;
	std::vector<double> mass;
	double i = 0;
	for (Star* star : stars) {
		index.push_back(i);
		i++;
		mass.push_back(star->mass);
	}
	InOut::write(index, mass, "initialConditionsMassBulge" + std::to_string(totalMass) + ".dat");
}

void Test::potentialSurfaceDensityBulge(){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 100; R < 30000; R = R + 100) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityBulge(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityBulge.dat");
}

void Test::potentialSurfaceDensityDisk(){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 100; R < 30000; R = R + 100) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityDisk(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityDisk.dat");
}

void Test::initialConditionsSampleDisk(){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	InitialConditions initialConditions = InitialConditions(&potential);
	double gridResolution = 0.001;
	Vec3D position = Vec3D(5, 5, 0);
	Vec3D volumeElement = Vec3D(gridResolution, gridResolution, gridResolution);
	double massInCell = potential.massDisk(position, volumeElement);
	int starID = 0;
	std::vector<Star*> starsInCell = initialConditions.massDisk(massInCell,starID); //stars with mass
	initialConditions.sampleDiskPositions(starsInCell, position, volumeElement);
	initialConditions.sampleDiskVelocities(starsInCell);
}

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
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
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
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	std::vector<Vec3D> Output;
	double border = 30000;
	ProgressBar progressBar = ProgressBar(-border, border-dx);
	std::cout << "Generating disk mass distribution" << std::endl;
	for (double x = -border; x < border; x += dx) {
		for (double y = -border; y < border; y += dx) {
			double starMass = potential.massBulge(Vec3D(x, y, z), Vec3D(dx, dx, 100));
			Output.push_back(Vec3D(x, y, starMass));
		}
		progressBar.Update(x);
		progressBar.Print();
	}

	InOut::write(Output, path+"massDistributionBulge_z" + std::to_string((int)z) + ".dat");

}

void Test::massDistributionTimer(){
	double totalMass = 0;
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	for (double x = 0; x < 10000; x = x + 100) {
		Vec3D corner = Vec3D(x, 0, -50);
		Vec3D volumeElement = Vec3D(10, 10, 10);
		Vec3D center = Vec3D::center(corner, corner+volumeElement);

		std::cout << "Volume Element at (center): " << center.length() * 1e-3 << "kpc" << std::endl;
		std::cout << "Volume Element dx:" << volumeElement.x * 1e-3 << "kpc" << std::endl;
		//std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
		double starMass = potential.massDisk(corner, volumeElement)+potential.massBulge(corner, volumeElement);
		totalMass += starMass;

		std::cout << "starMass: " << starMass << std::endl;

		//std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
		//std::cout << "Mass inside Volume: " << starMass << std::endl;
		//std::cout << "Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
		//startTime = std::chrono::steady_clock::now();
		//std::vector<Star*> stars = InitialConditions::massDisk(1e4);
		//endTime = std::chrono::steady_clock::now();
		//time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
		//std::cout << "Sample mass Disk - Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
	}
	std::cout << "totalMass: " << totalMass << std::endl;
}

void Test::velocityDispersionBulge(){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	for (double r = 100; r < 6000; r += 100) {
		std::cout << "r: " << r << " | velocityDistribution: " << potential.velocityDispersionBulge(r) << std::endl;
	}
}

void Test::velocityBulge(){
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	std::cout << "BulgeVelocityTest: start" << std::endl;
	std::vector<double> longitude;
	std::vector<double> dispersion;
	std::vector<double> bValues{ -4, -6, -8 };
	for (double b : bValues) {
		for (double l = -10; l < 10; l += 0.2) {
			longitude.push_back(l);
			std::vector<Vec3D*> velocities;
			double distance = 8490;
			double x = (8500 - distance * cos(b * 0.0174533) * cos(l * 0.0174533));
			double y = distance * cos(b * 0.0174533) * sin(l * 0.0174533);
			double z = distance * sin(b * 0.0174533);
			double r = gsl_hypot3(x, y, z);
			double R = gsl_hypot(x, y);
			double disp = potential.velocityDispersionBulge(r);
			//std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
			dispersion.push_back(disp);
			//std::cout << "radial dist:" << Potential::radialVelocityDispersionBulge(R, z) << std::endl;
		}
		InOut::write(longitude, dispersion, path + "/bulgeDispersion" + std::to_string((int)b)+".dat");
		longitude.clear();
		dispersion.clear();
	}
	Plot plot = Plot(path, path, true);
	plot.plot("bulgeDispersion", { });
}

void Test::velocityBulgeR() {
	std::vector<double> radius;
	std::vector<double> dispersion;
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	for (double r = 10; r < 1000; r += 1) {
		radius.push_back(r);
		double disp = potential.velocityDispersionBulge(r);
		std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
		dispersion.push_back(disp);
		//std::cout << "radial dist:" << Potential::radialVelocityDispersionBulge(R, z) << std::endl;
	}
	InOut::write(radius, dispersion, "bulgeDispersion.dat");

}

void Test::initialConditionsSampleBulgeVelocity(){
	Parameters parameters = Parameters();
	parameters.setNStars(1);
	Potential potential = Potential(&parameters);
	InitialConditions initialConditions = InitialConditions(&potential);
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initStars(starID, parameters.getNStars());
	double delta = 100;
	std::vector<double> bValues{ -4, -6, -8 };
	std::vector<double> longitude;
	std::vector<double> averageVelocity;
	for (double b : bValues) {
		for (double l = -10; l < 10; l += 0.2) {
			longitude.push_back(l);
			double distance = 8490;
			double x = (8500 - distance * cos(b * 0.0174533) * cos(l * 0.0174533));
			double y = distance * cos(b * 0.0174533) * sin(l * 0.0174533);
			double z = distance * sin(b * 0.0174533);
			double r = gsl_hypot3(x, y, z);
			initialConditions.sampleBulgePositions(stars, Vec3D(r, 0, 0), Vec3D(delta, delta, delta));
			initialConditions.sampleBulgeVelocities(stars);
			std::vector<Vec3D*> velocities;
			for (Star* star : stars) {
				velocities.push_back(&star->velocity);
			}
			//std::cout << "r: " << r << " | average velocity: " << Analysis::average(velocities) << " | Dispersion: " << Analysis::dispersion(velocities) << std::endl;
			averageVelocity.push_back(Analysis::average(velocities));
		}
		InOut::write(longitude, averageVelocity, path + "bulgeMeanVelocity" + std::to_string((int)b) + ".dat");
		longitude.clear();
		averageVelocity.clear();
	}
	Plot plot = Plot(path, path, true);
	plot.plot("bulgeMeanVelocity", { });
}

void Test::escapeVelocity(){
	double delta = 50;
	Parameters parameters = Parameters();
	Potential potential = Potential(&parameters);
	for(double x = 0; x < 10000; x += delta) {
		std::cout << "r:" << x << " | velocity: " << potential.escapeVelocity(&Vec3D(x, 0, 0)) << std::endl;
	}
}

void Test::initialConditionsInitFieldStars(){
	std::cout << "Testing field star creation" << std::endl;
	Potential potential = Potential(this);
	InitialConditions initialConditions = InitialConditions(&potential);
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initFieldStars(starID,focus,viewPoint,distance,dx,angle); //0.00029088833
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
	InOut::write(stars, path+"fieldStars.dat");
}
