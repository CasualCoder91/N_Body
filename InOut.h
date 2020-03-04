#pragma once
#include <fstream>
#include <vector>
#include <string>

#include "Star.h"
#include "Node.h"

namespace InOut{
	void write(std::vector<Star*> stars,std::string filename);
	void WriteWithLabel(std::vector<Star*> stars, std::string filename);
	void WriteAll(std::vector<Star*> stars, std::string filename);
	void write(Node* tree);
	void WriteRecursively(std::ofstream* file_ptr, Node* node_ptr);
	void write(std::vector<double> x, std::vector<double> y, std::string filename);
};

