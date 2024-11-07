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

	void printLoadedData();
	void calcElementJacobians();
	void calcElementsStiffnessMatrixes();

public:
	

	//Pattern methods
	void Run();
};

