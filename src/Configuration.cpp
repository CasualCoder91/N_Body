#include "Configuration.h"

Configuration::Configuration(const std::string& file) {
    Load(file);
}

void Configuration::Clear(){
    data.clear();
}

bool Configuration::Load(const std::string& file){
    std::ifstream inFile(file.c_str());

    if (!inFile.good()){
        std::cout << "Cannot read configuration file " << file << std::endl;
        return false;
    }

    while (inFile.good() && !inFile.eof())
    {
        std::string line;
        getline(inFile, line);

        // filter out comments
        if (!line.empty()){
            size_t pos = line.find('#');

            if (pos != std::string::npos){
                line = line.substr(0, pos);
            }
        }

        // split line into key and value
        if (!line.empty()){
            size_t pos = line.find('=');

            if (pos != std::string::npos){
                std::string key = Trim(line.substr(0, pos));
                std::string value = Trim(line.substr(pos + 1));

                if (!key.empty() && !value.empty()){
                    data[key] = value;
                }
            }
        }
    }

    return true;
}

bool Configuration::Contains(const std::string& key) const{
    return data.find(key) != data.end();
}

bool Configuration::Get(const std::string& key, std::string& value) const
{
    std::map<std::string, std::string>::const_iterator iter = data.find(key);

    if (iter != data.end()){
        value = iter->second;
        return true;
    }
    else{
        return false;
    }
}

bool Configuration::Get(const std::string& key, int& value) const{
    std::string str;

    if (Get(key, str)){
        value = atoi(str.c_str());
        return true;
    }
    else{
        return false;
    }
}

bool Configuration::Get(const std::string& key, long& value) const{
    std::string str;

    if (Get(key, str)){
        value = atol(str.c_str());
        return true;
    }
    else{
        return false;
    }
}

bool Configuration::Get(const std::string& key, double& value) const
{
    std::string str;

    if (Get(key, str)){
        value = atof(str.c_str());
        return true;
    }
    else{
        return false;
    }
}

double Configuration::GetDouble(const std::string& key) const
{
    std::string str;

    if (Get(key, str)) {
        return atof(str.c_str());
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        std::cout << "Using default value: 1" << std::endl;
        return 1;
    }
}

int Configuration::GetInt(const std::string& key) const{
    std::string str;

    if (Get(key, str)) {
        return atoi(str.c_str());
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        std::cout << "Using default value: 1" << std::endl;
        return 1;
    }
}

Vec3D Configuration::GetVec3D(const std::string& key) const{
    std::string str;
    Vec3D returnValue = Vec3D(0,0,0);
    if (Get(key, str)) {
        str = str.substr(1, str.size() - 2);
        std::stringstream ss(str);
        std::string item;
        std::getline(ss, item, ',');
        returnValue.x = std::stod(item);
        std::getline(ss, item, ',');
        returnValue.y = std::stod(item);
        std::getline(ss, item, ',');
        returnValue.z = std::stod(item);
        return returnValue;
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        std::cout << "Using default value: (0,0,0)" << std::endl;
        return returnValue;
    }
}

std::string Configuration::GetString(const std::string& key) const{
    std::map<std::string, std::string>::const_iterator iter = data.find(key);
    if (iter != data.end()) {
        return iter->second;
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        std::cout << "Using default value: NotFound" << std::endl;
        return "NotFound";
    }
}

bool Configuration::GetBool(const std::string& key) const{
    std::string str;

    if (Get(key, str)) {
        return (str == "true");
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        std::cout << "Using default value: false" << std::endl;
        return false;
    }
}

std::vector<double> Configuration::GetDoubleVector(const std::string& key) const{
    std::string str;
    std::vector<double> returnValue = {};
    if (Get(key, str)) {
        std::string temp;

        size_t i = 0, start = 0, end;
        do {
            end = str.find_first_of(',', start);
            temp = str.substr(start, end);
            if (isdigit(temp[0])){
                returnValue.push_back(atof(temp.c_str()));
                ++i;
            }
            start = end + 1;
        } while (start);
        return returnValue;
    }
    else {
        std::cout << key << " missing in cfg file!" << std::endl;
        returnValue.emplace_back(1);
        std::cout << "Using default value: {1}" << std::endl;
        return returnValue;
    }
}

bool Configuration::Get(const std::string& key, bool& value) const{
    std::string str;

    if (Get(key, str)){
        value = (str == "true");
        return true;
    }
    else{
        return false;
    }
}

bool Configuration::Get(const std::string& key, Vec3D& value) const{
    std::string str;
    if (Get(key, str)) {
        str = str.substr(1, str.size() - 2);
        std::stringstream ss(str);
        std::string item;
        std::getline(ss, item, ',');
        value.x = std::stod(item);
        std::getline(ss, item, ',');
        value.y = std::stod(item);
        std::getline(ss, item, ',');
        value.z = std::stod(item);
        return true;
    }
    else {
        return false;
    }
}

std::string Configuration::Trim(const std::string& str){
    size_t first = str.find_first_not_of(" \t");

    if (first != std::string::npos){
        size_t last = str.find_last_not_of(" \t");

        return str.substr(first, last - first + 1);
    }
    else{
        return "";
    }
}
