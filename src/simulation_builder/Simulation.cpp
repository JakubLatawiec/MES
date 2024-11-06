#include "Simulation.h"

#include <iostream>

void Simulation::printLoadedData()
{
    std::cout << "SimulationTime " << GlobalData.SimulationTime << std::endl;
    std::cout << "SimulationStepTime " << GlobalData.SimulationStepTime << std::endl;
    std::cout << "Conductivity " << GlobalData.Conductivity << std::endl;
    std::cout << "Alfa " << GlobalData.Alfa << std::endl;
    std::cout << "Tot " << GlobalData.Tot << std::endl;
    std::cout << "InitialTemp " << GlobalData.InitialTemp << std::endl;
    std::cout << "Density " << GlobalData.Density << std::endl;
    std::cout << "SpecificHeat " << GlobalData.SpecificHeat << std::endl;
    std::cout << "Nodes number " << GlobalData.nN << std::endl;
    std::cout << "Elements number " << GlobalData.nE << std::endl;

    std::cout << "*Node" << std::endl;
    for (int i = 0; i < GlobalData.nN; ++i) {
        std::cout << i + 1 << ": " << Grid.Nodes[i].x << ", " << Grid.Nodes[i].y << std::endl;
    }

    std::cout << "*Element, type=DC2D4" << std::endl;
    for (int i = 0; i < GlobalData.nE; ++i) {
        std::cout << i + 1 << ": ";
        for (int j = 0; j < 4; ++j) {
            std::cout << Grid.Elements[i].NodesID[j];
            if (j < 3) std::cout << ", ";
        }
        std::cout << std::endl;
    }
}

void Simulation::Run()
{
	std::cout << "RUN\n";
    printLoadedData();
}
