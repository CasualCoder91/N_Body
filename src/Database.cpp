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
		"viewPointZ REAL NOT NULL,"
		"bMcLuster INTEGER NOT NULL"
		");";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS star("
		"id INTEGER NOT NULL,"
		"id_simulation INTEGER NOT NULL,"
		"mass REAL,"
		"magnitude REAL,"
		"isCluster INTEGER,"
		"isObserved INTEGER NOT NULL,"
		"fkStar INTEGER,"
		"idCluster INTEGER, "
		"PRIMARY KEY (id),"
		"FOREIGN KEY (id_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS velocity("
		"id INTEGER PRIMARY KEY,"
		"x REAL,"
		"y REAL,"
		"z REAL,"
		"rH REAL,"
		"aHTP REAL,"
		"dHTP REAL,"
		"timestep INTEGER NOT NULL,"
		"id_star INTEGER NOT NULL)";
	this->exec(sql);
	sql = "CREATE UNIQUE INDEX IF NOT EXISTS idx_velocity_timestep_id_star "
		"ON velocity(timestep,id_star); ";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS position("
		"id INTEGER PRIMARY KEY,"
		"id_star INTEGER," //can be null because id_star of observed position initially unknown
		"x REAL,"
		"y REAL,"
		"z REAL,"
		"rH REAL,"
		"aHTP REAL,"
		"dHTP REAL,"
		"timestep INTEGER NOT NULL);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS timeStepAnalysis("
		"dt INTEGER NOT NULL,"
		"avgVel3DCluster REAL,"
		"avgVel2DCluster REAL,"
		"avgVel3DFS REAL,"
		"avgVel2DFS REAL,"
		"disp3DCluster REAL,"
		"disp2DCluster REAL,"
		"disp3DFS REAL,"
		"disp2DFS REAL,"
		"kinE REAL,"
		"potE REAL,"
		"totE REAL,"
		"id_simulation INTEGER NOT NULL,"
		"PRIMARY KEY(dt,id_simulation),"
		"FOREIGN KEY (id_simulation) "
			"REFERENCES simulation(id) "
			"ON DELETE CASCADE "
			"ON UPDATE NO ACTION);";
	this->exec(sql);
	sql = "CREATE TABLE IF NOT EXISTS powerLaw("
		"id_simulation INTEGER NOT NULL,"
		"position INTEGER NOT NULL,"
		"massLimit REAL,"
		"exponent REAL,"
		"PRIMARY KEY(id_simulation,position),"
		"FOREIGN KEY (id_simulation) "
		"REFERENCES simulation(id) "
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

int Database::insertSimulation(){
	if (!this->isOpen)
		this->open();
	std::string sql = "INSERT INTO simulation (n_stars,boxLength,dt,n_timesteps,title,outputTimestep,softening,precission,"
		"offsetX,offsetY,offsetZ,angle,distance,focusX,focusY,focusZ,viewPointX,viewPointY,viewPointZ,bMcLuster)"
		"VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11,?12,?13,?14,?15,?16,?17,?18,?19,?20)";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, Constants::nStars);
	sqlite3_bind_double(st, 2, Constants::boxLength);
	sqlite3_bind_double(st, 3, Constants::dt);
	sqlite3_bind_int(st, 4, Constants::nTimesteps);
	char* cstr = new char[Constants::title.length() + 1];
	strcpy_s(cstr, sizeof cstr, Constants::title.c_str());
	sqlite3_bind_text(st, 5, cstr, -1, SQLITE_TRANSIENT);
	sqlite3_bind_int(st, 6, Constants::outputTimestep);
	sqlite3_bind_double(st, 7, Constants::softening);
	sqlite3_bind_double(st, 8, Constants::precission);
	sqlite3_bind_double(st, 9, Constants::clusterLocation.x);
	sqlite3_bind_double(st, 10, Constants::clusterLocation.y);
	sqlite3_bind_double(st, 11, Constants::clusterLocation.z);
	sqlite3_bind_double(st, 12, Constants::angleOfView);
	sqlite3_bind_double(st, 13, Constants::distance);
	sqlite3_bind_double(st, 14, Constants::focus.x);
	sqlite3_bind_double(st, 15, Constants::focus.y);
	sqlite3_bind_double(st, 16, Constants::focus.z);
	sqlite3_bind_double(st, 17, Constants::viewPoint.x);
	sqlite3_bind_double(st, 18, Constants::viewPoint.y);
	sqlite3_bind_double(st, 19, Constants::viewPoint.z);
	sqlite3_bind_int(st, 20, Constants::bMcLuster);

	int returnCode = sqlite3_step(st);
	if (returnCode != SQLITE_DONE){
		throw "Could not insert new simulation";
	}
	sqlite3_finalize(st);
	int simulationID = getLastID();

	if (!Constants::bMcLuster) {
		this->insertPowerLaw(simulationID, Constants::massLimits, Constants::exponents);
	}

	return simulationID;
}

void Database::insertStars(int simulationID, std::vector<Star*>& stars, int timestep, bool clusterStars){
	std::cout << "Adding " << stars.size() << " stars to database" << std::endl;
	if (!this->isOpen)
		this->open();
	char* errorMessage;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string buffer = "INSERT INTO star (id,mass,id_simulation,isCluster,isObserved) VALUES (?1,?2,?3,?4,0)";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, buffer.c_str(), static_cast<int>(buffer.size()), &stmt, NULL);
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

int Database::insertAnalysis(int simulationID){
	std::string sql = "INSERT OR REPLACE INTO analysis (id_simulation) VALUES " 
		"((select id_simulation from analysis where id_simulation = ?1))";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, simulationID);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
	return getLastID();
}

void Database::insertAnalysisdtEnergy(int analysisID, int dt, double kinE, double potE){
	std::string sql = "INSERT OR REPLACE INTO timeStepAnalysis (dt,kinE,potE,totE,id_simulation,"
		"avgVel3DCluster, avgVel2DCluster, avgVel3DFS, avgVel2DFS, disp3DCluster, disp2DCluster, disp3DFS, disp2DFS) VALUES "
		"(?1,?2,?3,?4,?5,"
		"(SELECT avgVel3DCluster FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT avgVel2DCluster FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT avgVel3DFS FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT avgVel2DFS FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT disp3DCluster FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT disp2DCluster FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT disp3DFS FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1),"
		"(SELECT disp2DFS FROM timeStepAnalysis WHERE id_simulation = ?5 AND dt = ?1))";
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

void Database::insertAnalysisdtVelocity3D(int analysisID, int dt, double avgVel3DCluster, double disp3DCluster, double avgVel3DFS, double disp3DFS){
	std::string sql = "INSERT OR REPLACE INTO timeStepAnalysis "
		              "(dt,avgVel3DCluster,disp3DCluster,avgVel3DFS,disp3DFS,id_simulation,"
		              "kinE,potE,totE,avgVel2DCluster,avgVel2DFS,disp2DCluster,disp2DFS) VALUES "
		              "(?1,?2,?3,?4,?5,?6,"
		              "(SELECT kinE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT potE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT totE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT avgVel2DCluster FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT avgVel2DFS FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT disp2DCluster FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
					  "(SELECT disp2DFS FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1))";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, dt);
	sqlite3_bind_double(st, 2, avgVel3DCluster);
	sqlite3_bind_double(st, 3, disp3DCluster);
	sqlite3_bind_double(st, 4, avgVel3DFS);
	sqlite3_bind_double(st, 5, disp3DFS);
	sqlite3_bind_int(st, 6, analysisID);
	int returnCode = sqlite3_step(st);
	sqlite3_finalize(st);
}

void Database::insertAnalysisdtVelocity2D(int analysisID, int dt, double avgVel2DCluster, double disp2DCluster, double avgVel2DFS, double disp2DFS){
	std::string sql = "INSERT OR REPLACE INTO timeStepAnalysis "
		"(dt,avgVel2DCluster,disp2DCluster,avgVel2DFS,disp2DFS,id_simulation,"
		"kinE,potE,totE,avgVel3DCluster,avgVel3DFS,disp3DCluster,disp3DFS) VALUES "
		"(?1,?2,?3,?4,?5,?6,"
		"(SELECT kinE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT potE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT totE FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT avgVel3DCluster FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT avgVel3DFS FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT disp3DCluster FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1),"
		"(SELECT disp3DFS FROM timeStepAnalysis WHERE id_simulation = ?6 AND dt = ?1))";
	sqlite3_stmt* st;
	sqlite3_prepare(db, sql.c_str(), -1, &st, NULL);
	sqlite3_bind_int(st, 1, dt);
	sqlite3_bind_double(st, 2, avgVel2DCluster);
	sqlite3_bind_double(st, 3, disp2DCluster);
	sqlite3_bind_double(st, 4, avgVel2DFS);
	sqlite3_bind_double(st, 5, disp2DFS);
	sqlite3_bind_int(st, 6, analysisID);
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
	sqlite3_prepare_v2(db, buffer, 69, &stmt, NULL);
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
	sqlite3_prepare_v2(db, buffer2, 69, &stmt, NULL);
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
}

void Database::generateHEQ(int simulationID){
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT position.id,position.x,position.y,position.z,velocity.id,velocity.x,velocity.y,velocity.z,"
		"simulation.angle,simulation.focusX,simulation.focusY,simulation.focusZ "
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"INNER JOIN simulation on simulation.id = star.id_simulation "
		"WHERE position.timestep = velocity.timestep AND star.id_simulation = ?1";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);

	//local variables
	double angle = 0;
	Vec3D focus, focusHEQ, pGCA, vGCA;

	//store rows for velocity2D table for better performance
	struct rowInsert {
		int fk;
		double x, y, z, aHTP, dHTP;
	};
	std::vector<rowInsert> velocitiesInsert;
	std::vector<rowInsert> positionsInsert;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Vec3D pLSR, vLSR, pHCA, vHCA, pHEQ, vHEQ;
		if (angle == 0) { // Only need to init/read these variables once
			angle = sqlite3_column_double(stmt, 8);
			focus = Vec3D(sqlite3_column_double(stmt, 9), sqlite3_column_double(stmt, 10), sqlite3_column_double(stmt, 11));
			Projection::GCAtoLSR(focus, vGCA, pLSR, vLSR);
			Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);
			Projection::HCAtoHEQ(pHCA, vHCA, focusHEQ, vHEQ);
		}
		pGCA = Vec3D(sqlite3_column_double(stmt, 1), sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3));
		vGCA = Vec3D(sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6), sqlite3_column_double(stmt, 7));
		Projection::GCAtoLSR(pGCA, vGCA, pLSR, vLSR);
		Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);
		Projection::HCAtoHEQ(pHCA, vHCA, pHEQ, vHEQ);
		//Projection::project(pGCA, vGCA, lookAt, viewPoint);

		//position = Projection::projectPosition(position, lookAt, viewPoint,angle);
		rowInsert positionInsert = { sqlite3_column_int(stmt, 0), pHEQ.x, pHEQ.y, pHEQ.z, pHEQ.y- focusHEQ.y, pHEQ.z-focusHEQ.z };
		positionsInsert.emplace_back(positionInsert);
		//velocity = Projection::projectVelocity(velocity, lookAt, angle);
		rowInsert velocityInsert = { sqlite3_column_int(stmt, 4), vHEQ.x, vHEQ.y, vHEQ.z,0,0};
		velocitiesInsert.emplace_back(velocityInsert);
	}
	sqlite3_finalize(stmt);

	char* errorMessage;
	//insert velocities
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string buffer = "UPDATE velocity SET rH=?2, aHEQ=?3, dHEQ=?4 WHERE id=?1";
	sqlite3_prepare_v2(db, buffer.c_str(), static_cast<int>(buffer.size()), &stmt, nullptr);
	for (rowInsert row : velocitiesInsert) {
		sqlite3_bind_int(stmt, 1, row.fk);
		sqlite3_bind_double(stmt, 2, row.x);
		sqlite3_bind_double(stmt, 3, row.y);
		sqlite3_bind_double(stmt, 4, row.z);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
			printf(errorMessage);
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

	//insert positions
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	buffer = "UPDATE position SET rH=?2, aHEQ=?3, dHEQ=?4, aHTP=?5, dHTP=?6 WHERE id=?1";
	sqlite3_prepare_v2(db, buffer.c_str(), static_cast<int>(buffer.size()), &stmt, NULL);
	for (rowInsert row : positionsInsert) {
		sqlite3_bind_int(stmt, 1, row.fk);
		sqlite3_bind_double(stmt, 2, row.x);
		sqlite3_bind_double(stmt, 3, row.y);
		sqlite3_bind_double(stmt, 4, row.z);
		sqlite3_bind_double(stmt, 5, row.aHTP);
		sqlite3_bind_double(stmt, 6, row.dHTP);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

}

void Database::generateHTP(int simulationID)
{
	if (!this->isOpen)
		this->open();

	std::string query = "SELECT position.id,position.x,position.y,position.z,velocity.id,velocity.x,velocity.y,velocity.z,"
		"simulation.focusX,simulation.focusY,simulation.focusZ "
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"INNER JOIN simulation on simulation.id = star.id_simulation "
		"WHERE position.timestep = velocity.timestep AND star.id_simulation = ?1";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);

	//local variables
	bool initialized = false;
	Vec3D focus, focusHCA, pGCA, vGCA;
	//Matrix rotationM;

	struct rowInsert {
		int fk;
		double x, y, z;
	};
	std::vector<rowInsert> positionsInsert;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		Vec3D pLSR, vLSR, pHCA, vHCA, pHTP, vHTP;
		if (!initialized) { // Only need to init/read these variables once
			focus = Vec3D(sqlite3_column_double(stmt, 8), sqlite3_column_double(stmt, 9), sqlite3_column_double(stmt, 10));
			Projection::GCAtoLSR(focus, vGCA, pLSR, vLSR);
			Projection::LSRtoHCA(pLSR, vLSR, focusHCA, vHCA);
			focusHCA = focusHCA.normalize();
			/*rotationM = Matrix::rotation(focusHCA, Vec3D(1, 0, 0));*/
			initialized = true;

		}
		pGCA = Vec3D(sqlite3_column_double(stmt, 1), sqlite3_column_double(stmt, 2), sqlite3_column_double(stmt, 3));
		vGCA = Vec3D(sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6), sqlite3_column_double(stmt, 7));
		Projection::GCAtoLSR(pGCA, vGCA, pLSR, vLSR);
		Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);

		//pHCA = rotationM * pHCA; //rotate Position such that x points towards focus
		Projection::HCAtoHTP(pHCA, pHTP, focusHCA);

		rowInsert positionInsert = { sqlite3_column_int(stmt, 0), pHTP.x, pHTP.y, pHTP.z};
		positionsInsert.emplace_back(positionInsert);

	}
	sqlite3_finalize(stmt);

	char* errorMessage;

	//insert positions
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string buffer = "UPDATE position SET rH=?2, aHTP=?3, dHTP=?4 WHERE id=?1";
	sqlite3_prepare_v2(db, buffer.c_str(), static_cast<int>(buffer.size()), &stmt, NULL);
	for (rowInsert row : positionsInsert) {
		sqlite3_bind_int(stmt, 1, row.fk);
		sqlite3_bind_double(stmt, 2, row.x);
		sqlite3_bind_double(stmt, 3, row.y);
		sqlite3_bind_double(stmt, 4, row.z);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("Commit Failed!\n");
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);
}

void Database::generateMagnitude(int simulationID){
	if (!this->isOpen)
		this->open();

	std::string query = "SELECT star.id, star.mass, position.rH "
		"FROM star INNER JOIN position on position.id_star = star.id "
		"INNER JOIN simulation on simulation.id = star.id_simulation "
		"WHERE position.timestep = 0 AND star.id_simulation = ?1";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);

	//store rows for velocity2D table for better performance
	struct rowInsert {
		int idStar;
		double magnitude;
	};
	std::vector<rowInsert> stars;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		int id = sqlite3_column_int(stmt, 0);
		double mass = sqlite3_column_double(stmt, 1);
		double distance = sqlite3_column_double(stmt, 2);

		double lum = luminosity(mass);
		double magnitude = apparentMagnitude(lum, distance);

		rowInsert starsInsert = { id, magnitude };
		stars.emplace_back(starsInsert);

	}
	sqlite3_finalize(stmt);

	char* errorMessage;
	//insert magnitude
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string buffer = "UPDATE star SET magnitude=?2 WHERE id=?1";
	sqlite3_prepare_v2(db, buffer.c_str(), static_cast<int>(buffer.size()), &stmt, nullptr);
	for (rowInsert row : stars) {
		sqlite3_bind_int(stmt, 1, row.idStar);
		sqlite3_bind_double(stmt, 2, row.magnitude);
		if (sqlite3_step(stmt) != SQLITE_DONE)
		{
			printf("generateMagnitude: Commit Failed!\n");
			printf(errorMessage);
		}

		sqlite3_reset(stmt);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmt);

}

void Database::insertPowerLaw(int simulationID, std::vector<double> massLimits, std::vector<double> exponents){
	if (!this->isOpen)
		this->open();

	if (exponents.size() == massLimits.size() - 1)
		exponents.emplace_back(0);
	int position = this->selectLastID("powerLaw");
	position++;
	char* errorMessage;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	char buffer[] = "INSERT INTO powerLaw (id_simulation,position,massLimit,exponent) VALUES (?1,?2,?3,?4)";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, buffer, 85, &stmt, NULL);
	for (int i = 0; i < massLimits.size();++i) {
		sqlite3_bind_int(stmt, 1, simulationID);
		sqlite3_bind_int(stmt, 2, position+i);
		sqlite3_bind_double(stmt, 3, massLimits[i]);
		sqlite3_bind_double(stmt, 4, exponents[i]);
		if (sqlite3_step(stmt) != SQLITE_DONE){
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

void Database::selectConstants(int ID){
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT title, n_stars, boxlength, dt, n_timesteps, outputTimestep, softening, precission," 
		"offsetX, offsetY, offsetZ, angle, distance, focusX, focusY, focusZ, viewPointX, viewPointY, viewPointZ, bMcLuster "
		"FROM simulation "
        "Where ID = " + std::to_string(ID);
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), 79, &stmt, nullptr);

	if (sqlite3_step(stmt) == SQLITE_ROW) {
		Constants::title = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
		Constants::nStars = sqlite3_column_int(stmt, 1);
		Constants::boxLength = sqlite3_column_double(stmt, 2);
		Constants::dt = sqlite3_column_double(stmt, 3);
		Constants::nTimesteps = sqlite3_column_int(stmt, 4);
		Constants::outputTimestep = sqlite3_column_int(stmt, 5);
		Constants::softening = sqlite3_column_double(stmt, 6);
		Constants::precission =	sqlite3_column_double(stmt, 7);
		Constants::clusterLocation.x = sqlite3_column_double(stmt, 8);
		Constants::clusterLocation.y = sqlite3_column_double(stmt, 9);
		Constants::clusterLocation.z = sqlite3_column_double(stmt, 10);
		Constants::angleOfView = sqlite3_column_double(stmt, 11);
		Constants::distance = sqlite3_column_double(stmt, 12);
		Constants::focus.x = sqlite3_column_double(stmt, 13);
		Constants::focus.y = sqlite3_column_double(stmt, 14);
		Constants::focus.z = sqlite3_column_double(stmt, 15);
		Constants::viewPoint.x = sqlite3_column_double(stmt, 16);
		Constants::viewPoint.y = sqlite3_column_double(stmt, 17);
		Constants::viewPoint.z = sqlite3_column_double(stmt, 18);
		Constants::bMcLuster = sqlite3_column_int(stmt, 19);
	}
	sqlite3_finalize(stmt);
	return;
}

bool Database::printSimulations(){
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT id,title,n_stars,boxlength,dt,n_timesteps,outputTimestep FROM simulation";

	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	bool containsSimulations = false;
	std::string availableSimulations = "Available Simulations:\n";
	while(sqlite3_step(stmt) == SQLITE_ROW) {
		containsSimulations = true;
		availableSimulations += "ID: " + std::to_string(sqlite3_column_int(stmt, 0)) + " | Title: " + reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1)) + '\n';
	}
	sqlite3_finalize(stmt);
	if (!containsSimulations) {
		std::cout << "Database does not contain any simulations!" << std::endl;
	}
	else {
		std::cout << availableSimulations;
	}
	return containsSimulations;
}

std::vector<Vec3D> Database::selectVelocities3D(int simulationID, int timestep, bool fieldStars, bool clusterStars){
	std::vector<Vec3D> velocities = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT x,y,z FROM velocity INNER JOIN star on velocity.id_star=star.id "
		                "WHERE star.id_simulation = " + std::to_string(simulationID);
	if (timestep != -1) {
		query += " AND velocity.timestep = " + std::to_string(timestep);
	}
	if (!(fieldStars && clusterStars)) { // if not field and cluster stars
		if (fieldStars)
			query += " AND star.isCluster = 0 ";
		else if (clusterStars)
			query += " AND star.isCluster = 1 ";
		else
			std::cout << "selectVelocities3D parameters make no sense. Selecting cluster and field stars" << std::endl;
	}
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		velocities.emplace_back(sqlite3_column_double(stmt, 0), sqlite3_column_double(stmt, 1), sqlite3_column_double(stmt, 2));
	}
	sqlite3_finalize(stmt);
	return velocities;
}

std::vector<Vec2D> Database::selectVelocitiesHTP(int simulationID, int timestep, bool fieldStars, bool clusterStars, double minMagnitude){
	std::vector<Vec2D> velocities = {};
	if (!this->isOpen)
		this->open();

	std::string query = "SELECT velocity.aHTP ,velocity.dHTP "
		"FROM star "
		"INNER JOIN velocity on velocity.id_star = star.id "
		"WHERE star.id_simulation =  " +std::to_string(simulationID);
	if (minMagnitude != -1) {
		query += " AND star.magnitude < " + std::to_string(minMagnitude);
	}
	if (timestep != -1) {
		query += " AND velocity.timestep = " + std::to_string(timestep);
	}

	if (!(fieldStars && clusterStars)) { // if not field and cluster stars
		if (fieldStars)
			query += " AND star.isCluster = 0 ";
		else if (clusterStars)
			query += " AND star.isCluster = 1 ";
		else
			std::cout << "selectVelocitiesHTP parameters make no sense. Selecting cluster and field stars" << std::endl;
	}
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		velocities.emplace_back(sqlite3_column_double(stmt, 0), sqlite3_column_double(stmt, 1));
	}
	sqlite3_finalize(stmt);
	return velocities;
}

std::vector<int> Database::selectTimesteps(int simulationID){
	std::vector<int> timeSteps = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT DISTINCT timestep FROM velocity INNER JOIN star on velocity.id_star=star.id "
		                "WHERE star.id_simulation = " + std::to_string(simulationID);
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

std::vector<Star> Database::selectStars(int simulationID, int timeStep, bool observed){
	std::vector<Star> stars = {};
	if (!this->isOpen)
		this->open();
	std::string query = "SELECT star.id,mass,position.x,position.y,position.z,velocity.x,velocity.y,velocity.z " 
		"FROM star INNER JOIN velocity on velocity.id_star = star.id "
		"INNER JOIN position on position.id_star = star.id "
		"where position.timestep = ?1 AND velocity.timestep = ?1 AND star.id_simulation = ?2 "
		"and star.isObserved = ?3";
	sqlite3_stmt* stmt;
	if (sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr) != SQLITE_OK) {
		printf("Database::selectStars error\n");
	}
	sqlite3_bind_int(stmt, 1, timeStep);
	sqlite3_bind_int(stmt, 2, simulationID);
	sqlite3_bind_int(stmt, 3, observed);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		stars.push_back(Star(sqlite3_column_int(stmt, 0),sqlite3_column_double(stmt,1), 
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

void Database::outputStars2D(int simulationID, std::string filePath, int timestep)
{
	if (!this->isOpen)
		this->open();

	std::vector<int> timesteps = { timestep };
	if (timestep == -1) {
		timesteps = selectTimesteps(simulationID);
	}

	for (int timestep : timesteps) {

		std::string query = "SELECT mass, position.rH, position.aHTP, position.dHTP "
			"FROM star INNER JOIN position on position.id_star = star.id "
			"WHERE star.id_simulation = ?1 "
			"AND position.timestep = ?2";
		sqlite3_stmt* stmt;
		sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
		sqlite3_bind_int(stmt, 1, simulationID);
		sqlite3_bind_int(stmt, 2, timestep);
		std::ofstream file(filePath+std::to_string(timestep)+".dat");
		/* dump columns names into the file */
		for (int i = 0; i < 4; i++) {
			file << sqlite3_column_name(stmt, i) << ',';
		}
		file << "\n";

		while (sqlite3_step(stmt) == SQLITE_ROW) {
			for (int i = 0; i < 4; i++) {
				file << sqlite3_column_text(stmt, i) << ',';
			}
			file << "\n";
		}
		sqlite3_finalize(stmt);
		file.close();
	}
}



std::vector<std::vector<Point>> Database::selectPoints(int simulationID, int timeStep, int nTimeSteps, double minMagnitude, bool observed){
	std::vector<std::vector<Point>> points;
	std::vector<Point> timeStepPoints;
	if (nTimeSteps < 1) {
		std::cout << "Amount of timesteps must be at least 1" << std::endl;
		return points;
	}
	std::string query = "SELECT star.id, position.timestep, position.aHTP, position.dHTP, velocity.aHTP, velocity.dHTP, star.isCluster, star.magnitude "
		"FROM star "
		"LEFT JOIN position on position.id_star = star.id "
		"LEFT JOIN velocity on velocity.id_star = star.id AND position.timestep = velocity.timestep "
		"WHERE star.id_simulation = ?1 AND position.timestep IN (?2,?3) AND position.timestep = velocity.timestep "
		"AND star.isObserved = ?4";
		if (minMagnitude != -1) {
			query += " AND star.magnitude < ?5";
		}
	query += " order by position.timestep";
	sqlite3_stmt* stmt;
	sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.size()), &stmt, nullptr);
	sqlite3_bind_int(stmt, 1, simulationID);
	sqlite3_bind_int(stmt, 2, timeStep);
	sqlite3_bind_int(stmt, 3, timeStep + nTimeSteps -1);
	sqlite3_bind_int(stmt, 4, observed);
	if (minMagnitude != -1) {
		sqlite3_bind_double(stmt, 5, minMagnitude);
	}

	int currentTimeStep = timeStep;

	while (sqlite3_step(stmt) == SQLITE_ROW) {
		if (sqlite3_column_int(stmt, 1) != currentTimeStep) {
			points.emplace_back(timeStepPoints);
			timeStepPoints.clear();
			currentTimeStep = sqlite3_column_int(stmt, 1);
		}
		timeStepPoints.emplace_back(sqlite3_column_int(stmt, 0), //id
			sqlite3_column_double(stmt, 2),  //x
			sqlite3_column_double(stmt, 3),  //y
			sqlite3_column_double(stmt, 4),  //vx
			sqlite3_column_double(stmt, 5),  //vy
			sqlite3_column_int(stmt,6),      //clusterStar
			sqlite3_column_double(stmt, 7)); //magnitude
	}
	points.emplace_back(timeStepPoints);//gotta insert Points at last timestep
	sqlite3_finalize(stmt);

	return points;
}

void Database::updatePoints(std::vector<std::vector<Point>>& points)
{
	int timeStep = 0;

	char* errorMessage;
	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string queryVelocity = "update velocity set aHTP=?1, dHTP=?2 where timestep = ?3 and id_star = ?4";
	sqlite3_stmt* stmtVelocity;
	sqlite3_prepare_v2(db, queryVelocity.c_str(), static_cast<int>(queryVelocity.size()), &stmtVelocity, nullptr);
	for (Point& point : points[timeStep]) {
		sqlite3_bind_double(stmtVelocity, 1, point.velocity[0]);
		sqlite3_bind_double(stmtVelocity, 2, point.velocity[1]);
		sqlite3_bind_int(stmtVelocity, 3, timeStep);
		sqlite3_bind_int(stmtVelocity, 4, point.id);
		if (sqlite3_step(stmtVelocity) != SQLITE_DONE)
		{
			printf("updatePoints: Commit Failed!\n");
			printf(errorMessage);
		}
		sqlite3_reset(stmtVelocity);

	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmtVelocity);
	std::cout << "Database: velocity.aHTP and velocity.dHTP updated" << std::endl;


	sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
	std::string queryStar = "update star set idCluster=?1 where id=?2";
	sqlite3_stmt* stmtStar;
	sqlite3_prepare_v2(db, queryStar.c_str(), static_cast<int>(queryStar.size()), &stmtStar, nullptr);
	for (Point& point : points[timeStep]) {
		sqlite3_bind_int(stmtStar, 1, point.cluster);
		sqlite3_bind_int(stmtStar, 2, point.id);
		if (sqlite3_step(stmtStar) != SQLITE_DONE)
		{
			printf("updatePoints: Commit Failed!\n");
			printf(errorMessage);
		}
		sqlite3_reset(stmtStar);
	}
	sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);
	sqlite3_finalize(stmtStar);
	std::cout << "Database: idCluster updated" << std::endl;
}

