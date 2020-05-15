#include "LookupTable.h"

LookupTable::LookupTable(std::string filename, std::string delimiter) {
    this->delimiter = delimiter;
    this->fileName = filename;
    init();
}

bool LookupTable::isEmpty(){
    return map.empty();
}

double LookupTable::get(double key){
    std::map<double, double>::iterator low, prev;
    low = map.lower_bound(key); // lowest value above used key
    if (low == map.end()) {
        std::cout << "UhOh! No matching value found in LookupTable.." << std::endl;
        std::cin.clear();
        std::cin.get();
        return -1;
    }
    else if (low == map.begin()) {
        return low->second;
    }
    else {
        prev = std::prev(low); //next lower value 
        if ((key - prev->first) < (low->first - key)) //if previous value closer to key
            return prev->second;
        else
            return low->second;
    }
}

void LookupTable::setMap(std::vector<double> keys, std::vector<double> values){
    if (keys.size() != values.size()) {
        throw  "Vector size must be equal";
    }
    for (int i = 0; i < keys.size(); i++) {
        map.insert(std::make_pair(keys.at(i), values.at(i)));
    }
}

void LookupTable::makeFile(std::string pFileName, std::string header){
    if (pFileName.size() != 0)
        fileName = pFileName;
    std::ofstream file(filePath + fileName);
    if (header.size()>0)
        file << header << '\n';
    //no NOT parallel this one
    for (std::map<double,double>::const_iterator it = map.begin();it != map.end(); ++it){
        file << it->first << ", " << it->second << '\n';
    }
    file.close();
}

void LookupTable::init(){

    std::string line;
    std::ifstream file(filePath+fileName);

    double key, value;
    while (std::getline(file, line)) {
        std::string firstToken = line.substr(0, line.find(delimiter));
        if (checkIsDouble(firstToken)) {
            size_t pos = line.find(delimiter);
            key = std::stod(line.substr(0, pos));
            line.erase(0, pos + delimiter.length());
            value = std::stod(line, nullptr);
            map.insert(std::make_pair(key, value));
        }
    }

    if (map.size() == 0 && debug) {
        std::cout << "Map initialization failed. File not found or corrupted:" << std::endl;
        std::cout << "Path: " << filePath << " | Filename: " << fileName << std::endl;
        std::cin.get();
    }

}


bool LookupTable::checkIsDouble(std::string inputString) {
    char* end;
    double result = strtod(inputString.c_str(), &end);
    if (end == inputString.c_str() || *end != '\0')
        return false;
    return true;
}