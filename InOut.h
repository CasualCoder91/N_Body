#pragma once
#include <fstream>
#include <vector>
#include <string>

#include "Star.h"
#include "Node.h"

namespace InOut{
	void Write(std::vector<Star*> stars,std::string filename);
	void WriteWithLabel(std::vector<Star*> stars, std::string filename);
	void WriteAll(std::vector<Star*> stars, std::string filename);
	void Write(Node* tree);
	void WriteRecursively(std::ofstream* file_ptr, Node* node_ptr);
};

