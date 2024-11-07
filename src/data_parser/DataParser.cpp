#include "DataParser.h"

#include <sstream>
#include <fstream>

std::vector<std::string> DataParser::split(const std::string& s, char delimeter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimeter)) {
        tokens.push_back(token);
    }
    return tokens;
}

void DataParser::ParseData(const std::string& filename, GlobalData& globalData, Grid& grid)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key;

        ss >> key;
        if (key == "SimulationTime") {
            ss >> globalData.SimulationTime;
        }
        else if (key == "SimulationStepTime") {
            ss >> globalData.SimulationStepTime;
        }
        else if (key == "Conductivity") {
            ss >> globalData.Conductivity;
        }
        else if (key == "Alfa") {
            ss >> globalData.Alfa;
        }
        else if (key == "Tot") {
            ss >> globalData.Tot;
        }
        else if (key == "InitialTemp") {
            ss >> globalData.InitialTemp;
        }
        else if (key == "Density") {
            ss >> globalData.Density;
        }
        else if (key == "SpecificHeat") {
            ss >> globalData.SpecificHeat;
        }
        else if (key == "Nodes") {
            ss.ignore(7);
            ss >> globalData.nN;
            grid.Nodes.resize(globalData.nN);
        }
        else if (key == "Elements") {
            ss.ignore(7);
            ss >> globalData.nE;
            grid.Elements.resize(globalData.nE);
        }
        else if (key == "*Node") {
            for (int i = 0; i < globalData.nN; ++i) {
                std::getline(file, line);
                std::vector<std::string> tokens = split(line, ',');
                grid.Nodes[i].x = std::stod(tokens[1]);
                grid.Nodes[i].y = std::stod(tokens[2]);
            }
        }
        else if (key == "*Element,") {
            for (int i = 0; i < globalData.nE; ++i) {
                std::getline(file, line);
                std::vector<std::string> tokens = split(line, ',');
                for (int j = 0; j < 4; ++j) {
                    //grid.Elements[i].m_NodesID[j] = std::stoi(tokens[j + 1]);
                    grid.Elements[i].setNodesID(j, std::stoi(tokens[j + 1]));
                }
            }
        }
        else if (key == "*BC") {
            break;
        }
    }

    file.close();
}
