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
            int pos = line.find('#');

            if (pos != std::string::npos){
                line = line.substr(0, pos);
            }
        }

        // split line into key and value
        if (!line.empty()){
            int pos = line.find('=');

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
        return 1;
    }
}




bool Configuration::Get(const std::string& key, bool& value) const
{
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
    int first = str.find_first_not_of(" \t");

    if (first != std::string::npos){
        int last = str.find_last_not_of(" \t");

        return str.substr(first, last - first + 1);
    }
    else{
        return "";
    }
}