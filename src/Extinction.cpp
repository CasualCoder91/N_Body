#include "..\include\Extinction.h"

Extinction::Extinction()
{
	//load Extinction map from file
    std::string line;
    std::ifstream file("src/LookupTables/extinction.dat");

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(tokens));
        std::vector<double> d_tokens;
        for (std::string token : tokens) {
            d_tokens.push_back(std::stod(token));
        }
        map.push_back(d_tokens);
    }
}

void Extinction::set_extinction(Star& star)
{
    Vec3D pLSR, vLSR, pHCA, vHCA, pHEQ, vHEQ;
    Projection::GCAtoLSR(star.position, star.velocity, pLSR, vLSR);
    Projection::LSRtoHCA(pLSR, vLSR, pHCA, vHCA);
    Projection::HCAtoHGP(pHCA, vHCA, pHEQ, vHEQ);

    for (size_t i = 0; i < map.size();++i) //lines
    {
        if (abs(pHEQ.y - map[i][0] * Constants::degInRad) < 0.25 && abs(pHEQ.z - map[i][1] * Constants::degInRad) < 0.25) 
        {
            //std::cout << map[i][0] << " " << map[i][1] << " " << map[i][4] << std::endl;
            double total_extinction = 0;
            for (size_t j = 1; j < map[i][4]*4; j=j+4) 
            {
                total_extinction += map[i][4 + j + 1]; // wird letste mag mitgenommen?
                //std::cout << map[i][4 + j] << std::endl;
                if (pHEQ.x/1000 < map[i][4 + j]) 
                {
                    star.extinction = total_extinction;
                    return;
                }
            }
        }
    }
    //std::cout << "no excitction found" << std::endl;
    //std::cin.clear();
    //std::cin.get();

    return;
}
