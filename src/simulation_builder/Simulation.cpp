#include "Simulation.h"

#include <iostream>

void Simulation::printLoadedData()
{
    std::cout << "LOADED DATA:\n";
    std::cout << "SimulationTime " << m_GlobalData.SimulationTime << std::endl;
    std::cout << "SimulationStepTime " << m_GlobalData.SimulationStepTime << std::endl;
    std::cout << "Conductivity " << m_GlobalData.Conductivity << std::endl;
    std::cout << "Alfa " << m_GlobalData.Alfa << std::endl;
    std::cout << "Tot " << m_GlobalData.Tot << std::endl;
    std::cout << "InitialTemp " << m_GlobalData.InitialTemp << std::endl;
    std::cout << "Density " << m_GlobalData.Density << std::endl;
    std::cout << "SpecificHeat " << m_GlobalData.SpecificHeat << std::endl;
    std::cout << "Nodes number " << m_GlobalData.nN << std::endl;
    std::cout << "Elements number " << m_GlobalData.nE << std::endl;

    std::cout << "*Node" << std::endl;
    for (int i = 0; i < m_GlobalData.nN; ++i) {
        std::cout << "\t" << i + 1 << ": " << m_Grid.Nodes[i].x << ", " << m_Grid.Nodes[i].y << std::endl;
    }

    std::cout << "*Element, type=DC2D4" << std::endl;
    for (int i = 0; i < m_GlobalData.nE; ++i) {
        std::cout << i + 1 << ": ";
        for (int j = 0; j < 4; ++j) {
            std::cout << m_Grid.Elements[i].getNodesID()[j];
            if (j < 3) std::cout << ", ";
        }
        std::cout << "\n";
    }

    std::cout << std::endl;
}

void Simulation::printIPC()
{
    std::cout << "INTEGRATION POINTS COUNT: " << m_IPC << "\n";
    std::cout << std::endl;
}

void Simulation::printElements()
{
    for (size_t i = 0; i < m_GlobalData.nE; ++i)
    {
        std::cout << "ELEMENT " << i + 1 << ":\n";
        Element e = m_Grid.Elements[i];
        e.PrintJacobianMatrixes();
        e.PrintStiffnessMatrixes();
        e.PrintStiffnessMatrix();
        std::cout << std::endl;
    }
}

void Simulation::calcElementJacobians()
{
    for (auto& element : m_Grid.Elements)
    {
        auto& elementNodes = m_Grid.getElementNodes(element);
        element.CalcJacobians(elementNodes);
    }
        
}

void Simulation::calcElementsStiffnessMatrixes()
{
    for (auto& element : m_Grid.Elements)
        element.CalcStiffnessMatrixes(m_GlobalData.Conductivity);
}

void Simulation::Run()
{
    //Debug loaded data
    printLoadedData();

    //Debug integration points count
    printIPC();

    //Calc universal element derivatives
    Element::CalcElementUniv(m_IPC);

    //Debug universal element derivatives
    Element::PrintElementUniv();

    //Calc Matrixes for every element
    calcElementJacobians();
    calcElementsStiffnessMatrixes();

    //Debug every element
    printElements();
}
