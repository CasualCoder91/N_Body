#include "Database.h"

char* Database::dataBaseDataPath = "Output/Database/Default.db";

Database::Database(){
	this->isOpen = false;
	this->open();
}

bool Database::open(char* pDataBaseDataPath){
	if (std::strlen(pDataBaseDataPath) == 0) {
		pDataBaseDataPath = this->dataBaseDataPath;
	}
	else {
		this->dataBaseDataPath =  pDataBaseDataPath;
	}
	int rc = sqlite3_open(this->dataBaseDataPath, &db);
	if (rc){
		std::cerr << "Error opening SQLite3 database: " << sqlite3_errmsg(db) << std::endl << std::endl;
		sqlite3_close(db);
		this->isOpen = false;
		return true;
	}
	this->isOpen = true;
	return false;
}

bool Database::exec(char* sql){
	if (this->open()) //error during opening connection
		return true;
	char* zErrMsg = 0;
	int rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);
	if (rc != SQLITE_OK) {
		fprintf(stderr, "SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
	}
	//sqlite3_close(db);
	return false;
}

void Database::setup(){
	if (this->open()) //error during opening connection
		return;
	char* sql = "CREATE TABLE IF NOT EXISTS simulation("
		"id INTEGER PRIMARY KEY,"
		"title TEXT NOT NULL,"
		"n_stars INTEGER NOT NULL,"
		"boxlength REAL NOT NULL,"
		"dt REAL NOT NULL,"
		"n_timesteps INTEGER NOT NULL,"
		"outputTimestep INTEGER NOT NULL,"
		"softening REAL NOT NULL,"
		"precission REAL NOT NULL,"
		"minMass REAL NOT NULL,"
		"maxMass REAL NOT NULL,"
		"alpha REAL NOT NULL,"
		"offsetX REAL NOT NULL,"
		"offsetY REAL NOT NULL,"
		"offsetZ REAL NOT NULL,"
		"angle REAL NOT NULL,"
		"distance REAL NOT NULL,"
		"focusX REAL NOT NULL,"
		"focusY REAL NOT NULL,"
		"focusZ REAL NOT NULL,"
		"viewPointX REAL NOT NULL,"
		"viewPointY REAL NOT NULL,"
		"viewPointZ REAL NOT NULL"
		");";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS star("
		"id INTEGER NOT NULL,"
		"id_simulation INTEGER NOT NULL,"
		"mass REAL NOT NULL,"
		"isCluster INTEGER NOT NULL,"
		"PRIMARY KEY (id),"
		"FOREIGN KEY (id_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS velocity("
		"id INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"z REAL NOT NULL,"
		"timestep INTEGER NOT NULL,"
		"id_star INTEGER NOT NULL,"
		"FOREIGN KEY (id_star) "
			"REFERENCES star(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS velocity2D("
		"fk_velocity INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"FOREIGN KEY (fk_velocity) "
			"REFERENCES velocity(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS position("
		"id INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"z REAL NOT NULL,"
		"timestep INTEGER NOT NULL,"
		"id_star INTEGER NOT NULL,"
		"FOREIGN KEY (id_star) "
			"REFERENCES star(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS position2D("
		"fk_position INTEGER PRIMARY KEY,"
		"x REAL NOT NULL,"
		"y REAL NOT NULL,"
		"FOREIGN KEY (fk_position) "
		"REFERENCES position(id) "
		"ON DELETE CASCADE "
		"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS analysis("
		"id_simulation INTEGER PRIMARY KEY,"
		"doEnergy INTEGER NOT NULL,"
		"doVelocity INTEGER NOT NULL,"
		"doVelocity2D INTEGER NOT NULL,"
		"FOREIGN KEY (id_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS timeStepAnalysis("
		"dt INTEGER NOT NULL,"
		"averageVelocity REAL,"
		"kinE REAL,"
		"potE REAL,"
		"totE REAL,"
		"id_analysis INTEGER NOT NULL,"
		"PRIMARY KEY(dt,id_analysis),"
		"FOREIGN KEY (id_analysis) "
			"REFERENCES analysis(id_simulation) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
}

int Database::getLastID(){
	std::string sql = "SELECT Last_Insert_Rowid()";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	int returnCode = sqlite3_step(st);
	int ID = -1;
	if (returnCode == SQLITE_ROW) {
		ID = sqlite3_column_int(st, 0);
	}
	sqlite3_finalize(st);
	return ID;
}

int Database::selectLastID(std::string table){
	if (this->open()) //error during opening connection
		return -1;
	std::string sql = "select id from "+table+" order by id desc LIMIT 1";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	int returnCode = sqlite3_step(st);
	int ID = sqlite3_column_int(st, 0);
	sqlite3_finalize(st);
	return ID;
}

int Database::insert(Parameters* parameters){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO simulation (n_stars,boxLength,dt,n_timesteps,title,outputTimestep,softening,precission,minMass,maxMass,alpha,"
		"offsetX,offsetY,offsetZ,angle,distance,focusX,focusY,focusZ,viewPointX,viewPointY,viewPointZ)"
		"VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20,?21,?22)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, parameters->getNStars());
	sqlite3_bind_double(st, 2, parameters->getBoxLength());
	sqlite3_bind_double(st, 3, parameters->getdt());
	sqlite3_bind_int(st, 4, parameters->getNTimesteps());
	char* cstr = new char[parameters->getTitle().length() + 1];
	std::strcpy(cstr, parameters->getTitle().c_str());
	sqlite3_bind_text(st, 5, cstr, -1, SQLITE_TRANSIENT);
	sqlite3_bind_int(st, 6, parameters->getOutputTimestep());
	sqlite3_bind_double(st, 7, parameters->getSoftening());
	sqlite3_bind_double(st, 8, parameters->getPrecission());
	sqlite3_bind_double(st, 9, parameters->getMinMass());
	sqlite3_bind_double(st, 10, parameters->getMaxMass());
	sqlite3_bind_double(st, 11, parameters->getAlpha());
	sqlite3_bind_double(st, 12, parameters->getClusterLocation().x);
	sqlite3_bind_double(st, 13, parameters->getClusterLocation().y);
	sqlite3_bind_double(st, 14, parameters->getClusterLocation().z);
	sqlite3_bind_double(st, 15, parameters->getAngle());
	sqlite3_bind_double(st, 16, parameters->getDistance());
	sqlite3_bind_double(st, 17, parameters->getFocus().x);
	sqlite3_bind_double(st, 18, parameters->getFocus().y);
	sqlite3_bind_double(st, 19, parameters->getFocus().z);
	sqlite3_bind_double(st, 20, parameters->getViewPoint().x);
	sqlite3_bind_double(st, 21, parameters->getViewPoint().y);
	sqlite3_bind_double(st, 22, parameters->getViewPoint().z);

	int returnCode = sqlite3_step(st);
	if (returnCode != SQLITE_DONE){
		throw "Could not insert new simulation";
	}
	sqlite3_finalize(st);
	int simulationID = getLastID();
	return simulationID;
}

void Database::insertStars(int simulationID, std::vector<Star*>& stars, int timestep, bool clusterStars){
	std::cout << "Adding stars to database" << std::endl;
	if (!this->isOpen)
		this->open();
	char* errorMessage;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer[] = "INSERT INTO star (id,mass,id_simulation,isCluster) VALUES (?1,?2,?3,?4)";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, buffer, strlen(buffer), &stmt, NULL);
	//#pragma omp parallel for
	for (int i = 0; i < stars.size();++i) {
		sqlite3_bind_int(stmt, 1, stars[i]->id);
		sqlite3_bind_double(stmt, 2, stars[i]->mass);
		sqlite3_bind_int(stmt, 3, simulationID);
		sqlite3_bind_int(stmt, 4, int(clusterStars));
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}
		sqlite3_reset(stmt);
		//insertStar(simulationID, stars[i], timestep);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

	Database::timestep(timestep, stars);

}

int Database::insertAnalysis(int simulationID, Analysis analysis){
	std::string sql = "INSERT OR REPLACE INTO analysis (id_simulation,doEnergy,doVelocity,doVelocity2D) VALUES " 
		"((select id_simulation from analysis where id_simulation = ?1),?2,?3,?4)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, simulationID);
	sqlite3_bind_int(st, 2, analysis.getbEnergy());
	sqlite3_bind_int(st, 3, analysis.getbAverageVelocity());
	sqlite3_bind_int(st, 4, analysis.getbAverage2DVelocity());
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
	return getLastID();
}

void Database::insertAnalysisdtEnergy(int analysisID, int dt, double kinE, double potE){
	std::string sql = "INSERT OR REPLACE INTO timeStepAnalysis (dt,averageVelocity,kinE,potE,totE,id_analysis) VALUES "
		"(?1,"
		"(select averageVelocity from timeStepAnalysis where dt = ?1),"
		"?2,"
		"?3,"
		"?4,"
		"?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, dt);
	sqlite3_bind_double(st, 2, kinE);
	sqlite3_bind_double(st, 3, potE);
	sqlite3_bind_double(st, 4, kinE+potE);
	sqlite3_bind_int(st, 5, analysisID);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::insertAnalysisdtVelocity(int analysisID, int dt, double velocity){
	std::string sql = "INSERT OR REPLACE INTO timeStepAnalysis (dt,averageVelocity,kinE,potE,totE,id_analysis) VALUES "
		"(?1,"
		"?2,"
		"(select kinE from timeStepAnalysis where dt = ?1),"
		"(select potE from timeStepAnalysis where dt = ?1),"
		"(select totE from timeStepAnalysis where dt = ?1),"
		"?3)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, dt);
	sqlite3_bind_double(st, 2, velocity);
	sqlite3_bind_int(st, 3, analysisID);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::timestep(int timestep, std::vector<Star*>& stars){
	if (!this->isOpen)
		this->open();
	char* errorMessage;

	//Insert velocities
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer[] = "INSERT INTO velocity (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, buffer, strlen(buffer), &stmt, NULL);
	for (unsigned i = 0; i < stars.size(); i++){
		sqlite3_bind_double(stmt, 1, stars[i]->velocity.x);
		sqlite3_bind_double(stmt, 2, stars[i]->velocity.y);
		sqlite3_bind_double(stmt, 3, stars[i]->velocity.z);
		sqlite3_bind_int(stmt, 4, stars[i]->id);
		sqlite3_bind_int(stmt, 5, timestep);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

	//insert positions
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer2[] = "INSERT INTO position (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_prepare_v2(db, buffer2, strlen(buffer2), &stmt, NULL);
	for (unsigned i = 0; i < stars.size(); i++) {
		sqlite3_bind_double(stmt, 1, stars[i]->position.x);
		sqlite3_bind_double(stmt, 2, stars[i]->position.y);
		sqlite3_bind_double(stmt, 3, stars[i]->position.z);
		sqlite3_bind_int(stmt, 4, stars[i]->id);
		sqlite3_bind_int(stmt, 5, timestep);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}
		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

	//#pragma omp parallel for
	//for (int i = 0; i < stars.size(); ++i) {
	//	insertVelocity(stars[i]->id,stars[i]->velocity, timestep);
	//	insertPosition(stars[i]->id, stars[i]->position, timestep);
	//}
}

void Database::generate2D(int simulationID){
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT position.id,position.x,position.y,position.z,velocity.id,velocity.x,velocity.y,velocity.z,"
		"simulation.angle,simulation.focusX,simulation.focusY,simulation.focusZ,simulation.viewPointX,simulation.viewPointY,simulation.viewPointZ "
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"INNER JOIN simulation on simulation.id = star.id_simulation "
		"WHERE position.timestep = velocity.timestep AND star.id_simulation = ?1";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);

	//local variables
	double angle = 0;
	Vec3D focus, viewPoint, lookAt, position, velocity;

	//store rows for velocity2D table for better performance
	struct rowInsert {
		int fk;
		double x, y;
	};
	std::vector<rowInsert> velocitiesInsert;
	std::vector<rowInsert> positionsInsert;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		if (angle == 0) { // Only need to init/read these variables once
			angle = sqlite3_column_double(stmt, 8);
			focus = Vec3D(sqlite3_column_double(stmt, 9), sqlite3_column_double(stmt, 10), sqlite3_column_double(stmt, 11));
			viewPoint = Vec3D(sqlite3_column_double(stmt, 12), sqlite3_column_double(stmt, 13), sqlite3_column_double(stmt, 14));
			lookAt = (focus - viewPoint).normalize();
		}
		position = Vec3D(sqlite3_column_double(stmt, 1), sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3));
		velocity = Vec3D(sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6), sqlite3_column_double(stmt, 7));
		Projection::project(position, velocity, lookAt, viewPoint);

		//position = Projection::projectPosition(position, lookAt, viewPoint,angle);
		rowInsert positionInsert = { sqlite3_column_int(stmt, 0),position.x,position.y };
		positionsInsert.emplace_back(positionInsert);
		//velocity = Projection::projectVelocity(velocity, lookAt, angle);
		rowInsert velocityInsert = { sqlite3_column_int(stmt, 4),velocity.x,velocity.y };
		velocitiesInsert.emplace_back(velocityInsert);

	}
	sqlite3_finalize(stmt);

	char* errorMessage;
	//insert velocities
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer[] = "INSERT INTO velocity2D (fk_velocity,x,y) VALUES (?1,?2,?3)";
	sqlite3_prepare_v2(db, buffer, strlen(buffer), &stmt, NULL);
	for (rowInsert row : velocitiesInsert) {
		sqlite3_bind_int(stmt, 1, row.fk);
		sqlite3_bind_double(stmt, 2, row.x);
		sqlite3_bind_double(stmt, 3, row.y);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

	//insert positions
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer2[] = "INSERT INTO position2D (fk_position,x,y) VALUES (?1,?2,?3)";
	sqlite3_prepare_v2(db, buffer2, strlen(buffer2), &stmt, NULL);
	for (rowInsert row : positionsInsert) {
		sqlite3_bind_int(stmt, 1, row.fk);
		sqlite3_bind_double(stmt, 2, row.x);
		sqlite3_bind_double(stmt, 3, row.y);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

}

void Database::insertStar(int simulationID, Star* star, int& timestep, bool clusterStar){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO star (id,mass,id_simulation,isCluster) VALUES (?1,?2,?3,?4)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, star->id);
	sqlite3_bind_double(st, 2, star->mass);
	sqlite3_bind_int(st, 3, simulationID);
	sqlite3_bind_int(st, 4, clusterStar);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
	insertVelocity(star->id, star->velocity, timestep);
	insertPosition(star->id, star->position, timestep);
}

void Database::insertPosition(int& idStar, Vec3D& position, int& timestep){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO position (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_double(st, 1, position.x);
	sqlite3_bind_double(st, 2, position.y);
	sqlite3_bind_double(st, 3, position.z);
	sqlite3_bind_int(st, 4, idStar);
	sqlite3_bind_int(st, 5, timestep);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::insertVelocity(int& idStar, Vec3D& velocity, int& timestep){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO velocity (x,y,z,id_star,timestep) VALUES (?1,?2,?3,?4,?5)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_double(st, 1, velocity.x);
	sqlite3_bind_double(st, 2, velocity.y);
	sqlite3_bind_double(st, 3, velocity.z);
	sqlite3_bind_int(st, 4, idStar);
	sqlite3_bind_int(st, 5, timestep);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

std::vector<SimulationData> Database::selectSimulationData(int ID){
	std::vector<SimulationData> simulations = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT id,title,n_stars,boxlength,dt,n_timesteps,outputTimestep FROM simulation";
	if (ID != -1) {
		query += " Where ID = " + std::to_string(ID);
	}
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		simulations.push_back(SimulationData(sqlite3_column_int(stmt, 0), reinterpret_cast<const char*>(sqlite3_column_text(stmt,1)), sqlite3_column_int(stmt, 2),
			sqlite3_column_double(stmt, 3), sqlite3_column_double(stmt, 4), sqlite3_column_int(stmt, 5), sqlite3_column_int(stmt, 6)));
	}
	sqlite3_finalize(stmt);
	return simulations;
}

std::vector<Vec3D> Database::selectVelocities(int timestep){
	std::vector<Vec3D> velocities = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT x,y,z FROM velocity ";
	if (timestep != -1) {
		query += " WHERE timestep = " + std::to_string(timestep);
	}
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		velocities.push_back(Vec3D(sqlite3_column_int(stmt, 0), sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2)));
	}
	sqlite3_finalize(stmt);
	return velocities;
}

std::vector<int> Database::selectTimesteps(){
	std::vector<int> timeSteps = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT DISTINCT timestep FROM velocity";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		timeSteps.push_back(sqlite3_column_int(stmt, 0));
	}
	sqlite3_finalize(stmt);
	return timeSteps;
}

std::vector<Star*> Database::selectStars(int simulationID, int timeStep){
	std::vector<Star*> stars = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT star.id,mass,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z " 
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"where position.timestep = ?1 AND velocity.timestep = ?1 AND star.id_simulation = ?2";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr) != SQLITE_OK) {
		// Error reporting and handling
	}
	sqlite3_bind_int(stmt, 1, timeStep);
	sqlite3_bind_int(stmt, 2, simulationID);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		stars.push_back(new Star(sqlite3_column_int(stmt, 0),sqlite3_column_double(stmt,1), 
			sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3), sqlite3_column_double(stmt, 4),
			sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6), sqlite3_column_double(stmt, 7)));
	}
	sqlite3_finalize(stmt);
	return stars;
}

void Database::outputStars(int simulationID, std::string filePath, bool allStars, bool clusterStars, bool fieldStars) {
	if (!this->isOpen)
		this->open();
	if (allStars){
		std::string query = "SELECT star.id,mass,position.timestep,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z "
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"WHERE star.id_simulation = ?1 "
		"AND position.timestep = velocity.timestep"; //dont know why but this is needed
		sqlite3_stmt* stmt;
		sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
		sqlite3_bind_int(stmt, 1, simulationID);
		std::ofstream file(filePath);
		/* dump columns names into the file */
		for (int i = 0; i < 9; i++) {
			file << sqlite3_column_name(stmt, i) << ',';
		}
		file << "\n";

		while (sqlite3_step(stmt) == SQLITE_ROW) {
			for (int i = 0; i < 9; i++) {
				file << sqlite3_column_text(stmt, i) << ',';
			}
			file << "\n";
		}
		sqlite3_finalize(stmt);
		file.close();
	}
	if (clusterStars) {
		std::string query = "SELECT star.id,mass,position.timestep,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z "
			"FROM star INNER JOIN velocity on velocity.id_star = star.id "
			"INNER JOIN position on position.id_star = star.id "
			"WHERE star.id_simulation = ?1 "
			"AND position.timestep = velocity.timestep " //dont know why but this is needed
			"AND isCluster = 1";
		sqlite3_stmt* stmt;
		sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
		sqlite3_bind_int(stmt, 1, simulationID);
		std::ofstream file(filePath+"Cluster.dat");

		/* dump columns names into the file */
		for (int i = 0; i < 9; i++) {
			file << sqlite3_column_name(stmt, i) << ',';
		}
		file << "\n";

		while (sqlite3_step(stmt) == SQLITE_ROW) {
			for (int i = 0; i < 9; i++) {
				file << sqlite3_column_text(stmt, i) << ',';
			}
			file << "\n";
		}
		sqlite3_finalize(stmt);
		file.close();
	}
	if (fieldStars) {
		std::string query = "SELECT star.id,mass,position.timestep,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z "
			"FROM star INNER JOIN velocity on velocity.id_star = star.id "
			"INNER JOIN position on position.id_star = star.id "
			"WHERE star.id_simulation = ?1 "
			"AND position.timestep = velocity.timestep " //dont know why but this is needed
			"AND isCluster = 0";
		sqlite3_stmt* stmt;
		sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
		sqlite3_bind_int(stmt, 1, simulationID);
		std::ofstream file(filePath + "Field.dat");

		/* dump columns names into the file */
		for (int i = 0; i < 9; i++) {
			file << sqlite3_column_name(stmt, i) << ',';
		}
		file << "\n";

		while (sqlite3_step(stmt) == SQLITE_ROW) {
			for (int i = 0; i < 9; i++) {
				file << sqlite3_column_text(stmt, i) << ',';
			}
			file << "\n";
		}
		sqlite3_finalize(stmt);
		file.close();
	}
}

std::vector<std::vector<Point>> Database::selectPoints(int simulationID, int timeStep, int nTimeSteps){
	std::vector<std::vector<Point>> points;
	std::vector<Point> timeStepPoints;
	if (nTimeSteps < 1) {
		std::cout << "Amount of timesteps must be at least 1" << std::endl;
		return points;
	}
	std::string query = "SELECT star.id,position.timestep,position2D.x,position2D.y "
		"FROM star "
		"INNER JOIN position on position.id_star = star.id "
		"INNER JOIN position2D on position.id = position2D.fk_position "
		"WHERE star.id_simulation = ?1 AND position.timestep IN (?2,?3) order by position.timestep";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);
	sqlite3_bind_int(stmt, 2, timeStep);
	sqlite3_bind_int(stmt, 3, timeStep + nTimeSteps -1);

	int currentTimeStep = timeStep;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		if (sqlite3_column_int(stmt, 1) != currentTimeStep) {
			points.emplace_back(timeStepPoints);
			timeStepPoints.clear();
			currentTimeStep = sqlite3_column_int(stmt, 1);
		}
		timeStepPoints.emplace_back(sqlite3_column_int(stmt, 0), sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3));

	}
	sqlite3_finalize(stmt);

	return points;
}

