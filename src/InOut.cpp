#include "InOut.h"

void InOut::write(std::vector<Star*> stars, std::string filename) {
	std::ofstream file(InOut::outputDirectory + filename);
	for (Star* star : stars) {
		file << star->position.print() << '\n';
	}
	file.close();
}

void InOut::writeWithLabel(std::vector<Star*> stars, std::string filename){
	std::ofstream file(InOut::outputDirectory + filename);
	for (int i = 0; i < stars.size();++i) {
		file << stars.at(i)->position.x << ',' << stars.at(i)->position.y << ',' << stars.at(i)->position.z << ',' << i << '\n';
	}
	file.close();
}

void InOut::writeAll(std::vector<Star*> stars, std::string filename){
	std::ofstream file(InOut::outputDirectory + filename);
	#pragma omp parallel for
	for (int i = 0; i < stars.size(); ++i) {
		file << "Star: " << i << '\n';
		file << stars.at(i)->dump() << '\n';
	}
	file.close();
}

void InOut::write(Node* tree, std::string filename){
	if (!tree->isRoot()) {
		throw "Write function may only be called on root node";
	}
	else {
		std::ofstream file(InOut::outputDirectory + filename);
		writeRecursively(&file, tree);
		file.close();
	}
	return;
}



void InOut::writeRecursively(std::ofstream* file_ptr,Node* node_ptr) {
	for (Node* child : node_ptr->children) {
		if (child) {
			InOut::writeRecursively(file_ptr, child);
		}
	}
	*file_ptr << node_ptr->print();
	return;
}

std::vector<Vec3D> InOut::readVectors(std::string filename){
	std::vector<Vec3D> vectors;
	std::string line;
	std::ifstream sVectors(filename);
	if (sVectors.is_open()){
		while (std::getline(sVectors, line)){
			std::cout << line << '\n';
		}
		sVectors.close();
	}
	return vectors;
}

std::vector<std::vector<double>> InOut::readDoubleMatrix(std::string filname){
	std::string line;
	std::ifstream file(filname);
	std::vector<double> row;
	std::vector<std::vector<double>> matrix;
	std::string delimiter = ",";

	while (std::getline(file, line)){
		std::string firstToken = line.substr(0, line.find(delimiter));
		if (checkIsDouble(firstToken)) {
			size_t pos = 0;
			while ((pos = line.find(delimiter)) != std::string::npos) {
				row.push_back(std::stod(line.substr(0, pos), nullptr));
				line.erase(0, pos + delimiter.length());
			}
			row.push_back(std::stod(line, nullptr));
			matrix.push_back(row);
			row.clear();
		}
	}

	return matrix;
}

bool InOut::checkIsDouble(std::string inputString){
	char* end;
	double result = strtod(inputString.c_str(), &end);
	if (end == inputString.c_str() || *end != '\0') 
		return false;
	return true;
}

void InOut::write(std::vector<double> x, std::vector<double> y, std::string filename, std::string header){
	if (x.size() != y.size()) {
		throw  "Vector size must be equal";
	}
	std::ofstream file(filename);
	if (header.size() > 0)
		file << header << '\n';
	//no NOT parallel this one
	for (int i = 0; i < x.size(); ++i) {
		file << x.at(i) <<','<< y.at(i) << '\n';
	}
	file.close();
}

void InOut::write(std::vector<Vec3D> line, std::string filename) {
	std::ofstream file(InOut::outputDirectory + filename);
	for (int i = 0; i < line.size(); ++i) {
		file << line.at(i).print() << '\n';
	}
	file.close();
}
