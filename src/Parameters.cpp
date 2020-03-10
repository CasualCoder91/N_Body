#include "Parameters.h"


Parameters::Parameters(){
    this->filePath = "./simulation.cfg";
    config.Load(filePath);
    if (!config.Get("n_Stars", n_Stars)) {
        std::cout << "n_Stars missing in " + filePath << std::endl;
    }
    if (!config.Get("boxLength", boxLength)) {
        std::cout << "boxLength missing in " + filePath << std::endl;
    }
    if (!config.Get("dt", dt)) {
        std::cout << "dt missing in " + filePath << std::endl;
    }
    if (!config.Get("nTimesteps", nTimesteps)) {
        std::cout << "nTimesteps missing in " + filePath << std::endl;
    }
    if (!config.Get("title", title)) {
        std::cout << "title missing in " + filePath << std::endl;
    }
    this->softening = 0.16;
    this->G = 4.483e-3;
}

double Parameters::getG() {
    return G;
}
double Parameters::getSoftening() {
    return softening;
}
int Parameters::getN_Stars() {
    return n_Stars;
}
double Parameters::getBoxLength() {
    return boxLength;
}
double Parameters::getdt() {
    return dt;
}
int Parameters::getNTimesteps() {
    return nTimesteps;
}
std::string Parameters::getTitle() {
    return title;
}
