#include "..\include\Test.h"

void Test::samplePotentialOutput(int nStars) {

	std::vector<Star*> diskStars = InitialConditions::initStars(0, nStars);
	std::vector<Star*> bulgeStars = InitialConditions::initStars(0, nStars);

	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

	InitialConditions initialConditions = InitialConditions();
	initialConditions.sampleDiskPositions(diskStars, Vec3D(-40000, -40000, -20000), Vec3D(80000, 80000, 40000));
	initialConditions.sampleBulgePositions(bulgeStars, Vec3D(-6000, -6000, -6000), Vec3D(12000, 12000, 12000));

	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Time needed: " << time_span.count() << "seconds" << std::endl;
	//positionsDisk.push_back(Potential::sampleDisk(-40, 40, -40, 40, -20, 20));
	//positionsBulge.push_back(Potential::sampleBuldge(-6, 6, -6, 6, -6, 6));
	InOut::write(diskStars, "potentialDiskPositionsSample"+std::to_string(nStars)+".dat");
	InOut::write(bulgeStars, "potentialBulgePositionsSample" + std::to_string(nStars) + ".dat");
}

void Test::potentialCircularVelocityOutput() {
	Potential potential = Potential();

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

	InOut::write(positions, velocities, "../src/Test/potentialVelocity.dat");
	InOut::write(positions, disk, "../src/Test/potentialVelocityDisk.dat");
	InOut::write(positions, blackHole, "../src/Test/potentialVelocityBlackHole.dat");
	InOut::write(positions, buldge, "../src/Test/potentialVelocityBuldge.dat");
	InOut::write(positions, halo, "../src/Test/potentialVelocityHalo.dat");

	std::string filename = "E:/Master_Thesis/VS_Project/N_Body/src/Test/potentialVelocity.py";
	std::string command = "python ";
	command += filename;
	system(command.c_str());

}

void Test::testfrequencyDistribution() {
	Potential potential = Potential();
	std::vector<Vec3D> Output;

	double z = 1; //kpc
	double delta = 500;
	for (double x = -30000; x < 30000; x += delta) {
		std::cout << "x=" << x << std::endl;
		for (double y = -30000; y < 30000; y += delta) {
			double starMass = potential.frequencyDistribution(Vec3D(x, y, z), Vec3D(delta, delta, delta));
			Output.push_back(Vec3D(x, y, starMass));
		}
	}

	InOut::write(Output, "frequencyDistribution_z" + std::to_string(z) + ".dat");

}

void Test::initialConditionsMassSalpeterOutput(int nStars) {
	InitialConditions testInitialConditions = InitialConditions();
	testInitialConditions.setNStars(nStars);
	int starID = 0;
	std::vector<Star*> stars = testInitialConditions.initStars(starID);
	testInitialConditions.initialMassSalpeter(stars, 0.08, 100);
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
	InitialConditions initialConditions = InitialConditions();
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
	Potential potential = Potential();
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 100; R < 30000; R = R + 100) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityBulge(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityBulge.dat");
}

void Test::potentialSurfaceDensityDisk(){
	Potential potential = Potential();
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 100; R < 30000; R = R + 100) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityDisk(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityDisk.dat");
}

void Test::initialConditionsSampleDisk(){
	Potential potential = Potential();
	InitialConditions initialConditions = InitialConditions();
	double gridResolution = 0.001;
	Vec3D position = Vec3D(5, 5, 0);
	Vec3D volumeElement = Vec3D(gridResolution, gridResolution, gridResolution);
	double massInCell = potential.massDisk(position, volumeElement);
	int starID = 0;
	std::vector<Star*> starsInCell = initialConditions.massDisk(massInCell,starID); //stars with mass
	initialConditions.sampleDiskPositions(starsInCell, position, volumeElement);
	initialConditions.sampleDiskVelocities(starsInCell);
}

void Test::massDistributionDiskOutput(double z){
	Potential potential = Potential();
	std::vector<Vec3D> Output;

	for (double x = -30; x < 30; x += 0.5) {
		std::cout << "x=" << x << std::endl;
		for (double y = -30; y < 30; y += 0.5) {
			double starMass = potential.massDisk(Vec3D(x, y, z), Vec3D(0.1, 0.1, 0.1));
			Output.push_back(Vec3D(x, y, starMass));
		}
	}

	InOut::write(Output, "massDistributionDisk_z" + std::to_string(z) + ".dat");
}

void Test::massDistributionBulgeOutput(double z){
	Potential potential = Potential();
	std::vector<Vec3D> Output;

	for (double x = -30; x < 30; x += 0.5) {
		std::cout << "x=" << x << std::endl;
		for (double y = -30; y < 30; y += 0.5) {
			double starMass = potential.massBulge(Vec3D(x, y, z), Vec3D(0.1, 0.1, 0.1));
			Output.push_back(Vec3D(x, y, starMass));
		}
	}

	InOut::write(Output, "massDistributionBulge_z" + std::to_string((int)z) + ".dat");
}

void Test::massDistributionTimer(){
	double totalMass = 0;
	Potential potential = Potential();
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

void Test::velocityDistributionTestBulge(){
	Potential potential = Potential();
	for (double r = 100; r < 6000; r += 100) {
		std::cout << "r: " << r << " | velocityDistribution: " << potential.velocityDistributionBulge(r) << std::endl;
	}
}

void Test::velocityBulgeTest(){
	Potential potential = Potential();
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
			double z = distance * sin(b * 0.0174533); // b = -4°
			double r = gsl_hypot3(x, y, z);
			double R = gsl_hypot(x, y);
			double disp = potential.velocityDistributionBulge(r);
			//std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
			dispersion.push_back(disp);
			//std::cout << "radial dist:" << Potential::radialVelocityDispersionBulge(R, z) << std::endl;
		}
		InOut::write(longitude, dispersion, "../src/Test/bulgeDispersion" + std::to_string((int)b)+".dat");
		longitude.clear();
		dispersion.clear();
	}
	std::string filename = "E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeDispersion.py";
	std::string command = "python ";
	command += filename;
	system(command.c_str());

}

void Test::velocityBulgeRTest() {
	std::vector<double> radius;
	std::vector<double> dispersion;
	Potential potential = Potential();
	for (double r = 10; r < 1000; r += 1) {
		radius.push_back(r);
		double disp = potential.velocityDistributionBulge(r);
		std::cout << "r: " << r << " | vel dist: " << disp << std::endl;
		dispersion.push_back(disp);
		//std::cout << "radial dist:" << Potential::radialVelocityDispersionBulge(R, z) << std::endl;
	}
	InOut::write(radius, dispersion, "bulgeDispersion.dat");

}

void Test::initialConditionsSampleBulgeVelocity(){
	Parameters testParameters = Parameters();
	InitialConditions initialConditions = InitialConditions();
	initialConditions.setNStars(1);
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initStars(0, starID);
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
			double z = distance * sin(b * 0.0174533); // b = -4°
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
		InOut::write(longitude, averageVelocity, "../src/Test/bulgeMeanVelocity" + std::to_string((int)b) + ".dat");
		longitude.clear();
		averageVelocity.clear();
	}
	std::string filename = "E:/Master_Thesis/VS_Project/N_Body/src/Test/bulgeMeanVelocity.py";
	std::string command = "python ";
	command += filename;
	system(command.c_str());
}

void Test::escapeVelocityTest(){
	double delta = 50;
	Potential potential = Potential();
	for(double x = 0; x < 10000; x += delta) {
		std::cout << "r:" << x << " | velocity: " << potential.escapeVelocity(&Vec3D(x, 0, 0)) << std::endl;
	}
}

void Test::initialConditionsInitFieldStars(){
	std::cout << "Testing field star creation" << std::endl;
	InitialConditions initialConditions = InitialConditions();
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	int starID = 0;
	std::vector<Star*> stars = initialConditions.initFieldStars(starID); //0.00029088833
	std::chrono::steady_clock::time_point endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "Time needed for calucation: " << time_span.count() << "seconds" << std::endl;
	InOut::write(stars, "fieldStars.dat");
}
