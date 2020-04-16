#include "..\include\Test.h"

void Test::samplePotentialOutput(int nStars) {
	Potential potential = Potential(Vec3D(0, 0, 0));

	std::vector<Vec3D> positionsDisk;
	std::vector<Vec3D> positionsBulge;
	for (int i = 0; i < nStars; i++) {
		positionsDisk.push_back(potential.sampleDisk(-40, 40, -40, 40, -20, 20));
		positionsBulge.push_back(potential.sampleBuldge(-6, 6, -6, 6, -6, 6));
	}
	InOut::write(positionsDisk, "potentialDiskPositionsSample"+std::to_string(nStars)+".dat");
	InOut::write(positionsBulge, "potentialBulgePositionsSample" + std::to_string(nStars) + ".dat");
}

void Test::potentialVelocityOutput() {
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<double> positions;
	std::vector<double> velocities;
	std::vector<double> disk;
	std::vector<double> blackHole;
	std::vector<double> buldge;
	std::vector<double> halo;

	for (double i = 0.1; i < 30; i += 0.1) {
		positions.push_back(i);
		Vec3D position = Vec3D(i, 0, 0);
		velocities.push_back(potential.circularVelocity(&position));
		disk.push_back(potential.circularVelocityDisk(&position));
		blackHole.push_back(potential.circularVelocityBlackHole(&position));
		buldge.push_back(potential.circularVelocityBulge(&position));
		halo.push_back(potential.circularVelocityHalo(&position));
	}

	InOut::write(positions, velocities, "potentialTest.dat");
	InOut::write(positions, disk, "potentialTestDisk.dat");
	InOut::write(positions, blackHole, "potentialTestBlackHole.dat");
	InOut::write(positions, buldge, "potentialTestBuldge.dat");
	InOut::write(positions, halo, "potentialTestHalo.dat");
}

void Test::testfrequencyDistribution() {
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<Vec3D> Output;

	double z = 1; //kpc

	for (double x = -30; x < 30; x += 0.5) {
		std::cout << "x=" << x << std::endl;
		for (double y = -30; y < 30; y += 0.5) {
			double starMass = potential.frequencyDistribution(Vec3D(x, y, z), Vec3D(0.1, 0.1, 0.1));
			Output.push_back(Vec3D(x, y, starMass));
		}
	}

	InOut::write(Output, "frequencyDistribution_z" + std::to_string(z) + ".dat");

}

void Test::initialConditionsMassSalpeterOutput(int nStars) {
	Parameters testParameters = Parameters();
	InitialConditions testInitialConditions = InitialConditions(&testParameters);
	testInitialConditions.setNStars(nStars);
	std::vector<Star*> stars = testInitialConditions.initStars(0);
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
	Parameters parameters = Parameters();
	InitialConditions initialConditions = InitialConditions(&parameters);
	std::vector<Star*> stars = initialConditions.initialMassBulge(totalMass);
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
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 0.1; R < 30; R = R + 0.1) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityBulge(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityBulge.dat");
}

void Test::potentialSurfaceDensityDisk(){
	Potential potential = Potential(Vec3D(0, 0, 0));
	std::vector<double> surfaceDensity;
	std::vector<double> radius;
	for (double R = 0.1; R < 30; R = R + 0.1) {
		radius.push_back(R);
		surfaceDensity.push_back(potential.surfaceDensityDisk(R));
	}
	InOut::write(radius, surfaceDensity, "potentialSurfaceDensityDisk.dat");
}

void Test::massDistributionDiskOutput(double z){
	Potential potential = Potential(Vec3D(0, 0, 0));
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
	Potential potential = Potential(Vec3D(0, 0, 0));
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
