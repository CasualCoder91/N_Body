#include "Simulation.h"

Simulation::Simulation(int id){
	this->id = id;
	this->database = new Database(Constants::database_path);
	this->potential = new MWPotential();
}

Simulation::Simulation(int id, Database* database){
	this->id = id;
	this->potential = new MWPotential();
	this->database = database;
}

void Simulation::setID(int id){
	this->id = id;
}

int Simulation::getID(){
	return id;
}

void Simulation::run(bool reuse_cluster){

	bool euler = false;

	reuse_cluster = reuse_cluster && this->id != 1; // can only reuse the cluster if it already exists.

	//Init stars
	InitialConditions initialConditions = InitialConditions(Simulation::potential);

	//Init clusterStars
	int nextStarIndex = database->selectLastID("star") + 1;
	size_t cluster_star_startindex = nextStarIndex;
	std::vector<Star> clusterStars;
	if (reuse_cluster) {
		clusterStars = database->select_stars(1, 0, false, 1);
		for (Star& star : clusterStars) {
			star.id = star.id + cluster_star_startindex;
		}
	}
	else {
		if (Constants::bMcLuster) {
			clusterStars = InOut::readMcLuster(nextStarIndex, Constants::mcluster_filepath);
		}
		else {
			clusterStars = initialConditions.initStars(nextStarIndex, Constants::nStars);
			double totalMass = initialConditions.brokenPowerLaw(clusterStars, Constants::massLimits, Constants::exponents);
			initialConditions.plummerSphere(clusterStars, totalMass, Constants::plummer_radius, Constants::G);
		}
		initialConditions.offsetCluster(clusterStars, Constants::clusterLocation);
		double circVel = potential->circularVelocity(&Constants::clusterLocation);
		Vec3D clusterVelocity = Vec3D(0, -circVel, 0);
		//initialConditions.sampleDiskVelocity(clusterVelocity, clusterLocation);
		for (Star& star : clusterStars) {
			star.velocity += clusterVelocity;
		}
	}
	database->insertStars(this->getID(), clusterStars, 0, true);

	//Init field stars
	nextStarIndex = database->selectLastID("star") + 1;
	std::vector<Star> fieldStars = initialConditions.initFieldStars(nextStarIndex, Constants::focus, Constants::viewPoint, Constants::distance, Constants::angleOfView);
	database->insertStars(this->getID(), fieldStars, 0, false);

	std::cout << "Starting integration" << std::endl;
	//Integrate
	Integrator integrator = Integrator(Constants::dt);
	//std::cout << "dt = " << integrator.dt << std::endl;
	//std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	Vec3D tlf = Vec3D(), brb = Vec3D();
	ProgressBar progressBar = ProgressBar(0, Constants::nTimesteps);

	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point endTime;

	for (int i = 0; i <= Constants::nTimesteps; i++) {
		if (i>0 && i % Constants::outputTimestep == 0) {
			if (reuse_cluster) {
				clusterStars = database->select_stars(1, i / Constants::outputTimestep, false, 1);
				for (Star& star : clusterStars) {
					star.id = star.id + cluster_star_startindex;
				}
			}
			database->timestep(i/Constants::outputTimestep, clusterStars);
			database->timestep(i/Constants::outputTimestep, fieldStars);
			progressBar.Update(i);
			progressBar.Print();
			if (Constants::nTimesteps - i < Constants::outputTimestep)
				break;
		}

		if (clusterStars.size() > 0 && !reuse_cluster) {
			Node::findCorners(tlf, brb, clusterStars);
			Node root = Node(tlf, brb, nullptr); // Cleanup/Destructor of tree needed in every timestep

			for (Star& star : clusterStars) {
				root.insert(&star);
			}
			root.calculateMassDistribution();

			//Force clusterStars -- use if Euler
			if (euler) {
				for (int j = 0; j < clusterStars.size(); ++j) {
					clusterStars[j].acceleration.reset();
					potential->applyForce(&clusterStars[j]);
					root.applyForce(clusterStars[j].position, clusterStars[j].acceleration);
				}
				integrator.euler(clusterStars);
			}
			else {
				integrator.Leapfrog(clusterStars, &root, this->potential);
			}
		}

		//Force fieldStars -- use if Euler
		if (euler) {
			if (fieldStars.size() > 0) {
				for (int i = 0; i < fieldStars.size(); ++i) {
					fieldStars[i].acceleration.reset();
					potential->applyForce(&fieldStars[i]);
				}
			}
		}

		if (fieldStars.size() > 0) {
			if (euler) {
				integrator.euler(fieldStars);
			}
			else {
				integrator.Leapfrog(fieldStars, this->potential);
			}
		}

	}
	endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	//InOut::write(clusterStars,"stars.dat");
	//InOut::write(&root);
	std::cout << "Time needed for integration: " << time_span.count() << "seconds" << std::endl;
	std::cout << "done" << std::endl;
}

Simulation::~Simulation(){
	delete potential;
}
