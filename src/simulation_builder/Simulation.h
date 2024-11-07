#pragma once

#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

class Simulation
{
private:
	//Pattern variables
	Simulation() = default;
	friend class SimulationBuilder;

	//Simulation variables
	int m_IPC{};
	GlobalData m_GlobalData{};
	Grid m_Grid{};

	//Simulation calculation methods
	void calcElementJacobians();
	void calcElementsStiffnessMatrixes();

	//Debug methods
	void printLoadedData();
	void printIPC();
	void printElements();

public:
	//Pattern methods
	void Run();
};

