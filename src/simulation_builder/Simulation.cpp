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
        std::cout << "\t" << i + 1 << ": " << m_Grid.Nodes[i].x << ", " << m_Grid.Nodes[i].y << (m_Grid.Nodes[i].isBorderCondition ? " (BC)" : "") << std::endl;
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
        //e.PrintJacobianMatrixes();
        //e.PrintStiffnessMatrixes();
        //e.PrintHbcMatrixes();
        //e.PrintPVector();
        e.PrintStiffnessMatrix();
        e.PrintCMatrixes();
        e.PrintCMatrix();
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

void Simulation::calcElementHbcMatrixes()
{
    for (auto& element : m_Grid.Elements)
    {
        auto& elementNodes = m_Grid.getElementNodes(element);
        element.CalcHbcMatrixes(elementNodes, m_GlobalData.Alfa);
    }
}

void Simulation::calcElementPVector()
{
    for (auto& element : m_Grid.Elements)
    {
        auto& elementNodes = m_Grid.getElementNodes(element);
        element.CalcPVector(elementNodes, m_GlobalData.Alfa, m_GlobalData.Tot);
    }
}

void Simulation::calcGlobalStifnessMatrix()
{
    m_EquationsSolver = EquationsSolver(m_GlobalData.nN);
    for (auto& element : m_Grid.Elements)
    {
        auto& nodeIDs = element.getNodesID();
        auto& stiffnessMatrix = element.getStifnessMatrix();
        size_t stiffnessMatrixRow = 0;
        size_t stiffnessMatrixCol = 0;

        for (auto row_ID : nodeIDs)
        {
            for (auto col_ID : nodeIDs)
            {
                m_EquationsSolver.setGlobalStiffnessMatrix(row_ID - 1, col_ID - 1, stiffnessMatrix(stiffnessMatrixRow, stiffnessMatrixCol));
                ++stiffnessMatrixCol;
            }
            ++stiffnessMatrixRow;
            stiffnessMatrixCol = 0;
        }
    }
}

void Simulation::calcGlobalPVector()
{
    for (auto& element : m_Grid.Elements)
    {
        auto& nodeIDs = element.getNodesID();
        auto& pVector = element.getPVector();

        size_t pVectorIndex = 0;
        for (auto nodeID : nodeIDs)
        {
            m_EquationsSolver.setGlobalPVector(nodeID - 1, pVector(pVectorIndex, 0));
            ++pVectorIndex;
        }
    }
}

void Simulation::calcCMatrixes()
{
    for (auto& element : m_Grid.Elements)
        element.CalcCMatrix(m_GlobalData.SpecificHeat, m_GlobalData.Density);
}

void Simulation::calcGlobalCMatrix()
{
    for (auto& element : m_Grid.Elements)
    {
        auto& nodeIDs = element.getNodesID();
        auto& cMatrix = element.getCMatrix();
        size_t cMatrixRow = 0;
        size_t cMatrixCol = 0;

        for (auto row_ID : nodeIDs)
        {
            for (auto col_ID : nodeIDs)
            {
                m_EquationsSolver.setGlobalCMatrix(row_ID - 1, col_ID - 1, cMatrix(cMatrixRow, cMatrixCol));
                ++cMatrixCol;
            }
            ++cMatrixRow;
            cMatrixCol = 0;
        }
    }
}

void Simulation::Run()
{
    Gauss::Initialize(3);
    Element::Initialize(3, 4);

    //Debug loaded data
    printLoadedData();

    //Debug integration points count
    printIPC();

    //Calc Matrixes for every element
    calcElementJacobians();
    calcElementsStiffnessMatrixes();

    //Calc border conditions for every element
    calcElementHbcMatrixes();
    calcElementPVector();

    //Calc global stiffness matrix
    calcGlobalStifnessMatrix();

    //Calc global P vector
    calcGlobalPVector();

    //Calc C matrixes
    calcCMatrixes();

    calcGlobalCMatrix();

    //Debug every element
    printElements();

    //Print stiffness global matrix
    std::cout << "\nGLOBAL H MATRIX:\n";
    m_EquationsSolver.PrintGlobalStiffnessMatrix();

    //Print global P vector
    std::cout << "\nGLOBAL P VECTOR:\n";
    m_EquationsSolver.PrintGlobalPVector();

    //Print global C matrix
    std::cout << "\nGLOBAL C MATRIX:\n";
    m_EquationsSolver.PrintGlobalCMatrix();

    //Solve equation
    m_EquationsSolver.SolveEquation(m_GlobalData.SimulationStepTime, m_GlobalData.SimulationTime, m_GlobalData.InitialTemp);

    //Print global T vector
    std::cout << "\nGLOBAL T VECTOR:\n";
    m_EquationsSolver.PrintGlobalTVector();
}
