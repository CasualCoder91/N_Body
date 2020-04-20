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
	std::ifstream sVectors("filename");
	if (sVectors.is_open()){
		while (std::getline(sVectors, line)){
			std::cout << line << '\n';
		}
		sVectors.close();
	}
	return vectors;
}

void InOut::write(std::vector<double> x, std::vector<double> y, std::string filename){
	if (x.size() != y.size()) {
		throw  "Vector size must be equal";
	}
	std::ofstream file(InOut::outputDirectory+filename);
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
