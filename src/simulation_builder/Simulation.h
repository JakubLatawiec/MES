#pragma once

#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"
#include "../data_types/EquationsSolver.h"

class Simulation
{
private:
	//Pattern variables
	Simulation() = default;
	friend class SimulationBuilder;

	//Simulation variables
	int m_IPC{};
	int m_SurfaceIPC{};
	GlobalData m_GlobalData{};
	Grid m_Grid{};
	EquationsSolver m_EquationsSolver{};

	//Simulation calculation methods
	void calcElementJacobians();
	void calcElementsStiffnessMatrixes();
	void calcElementHbcMatrixes();
	void calcElementPVector();
	void calcGlobalStifnessMatrix();
	void calcGlobalPVector();
	void calcCMatrixes();
	void calcGlobalCMatrix();

	//Debug methods
	void printLoadedData();
	void printIPC();
	void printElements();

public:
	//Pattern methods
	void Run();
};

