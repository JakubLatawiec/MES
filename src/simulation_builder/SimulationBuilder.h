#pragma once

#include "../data_parser/DataParser.h"
#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

#include "Simulation.h"

class SimulationBuilder
{
private:
	//Pattern variables
	Simulation m_Simulation;

public:
	//Configuration setters
	SimulationBuilder& LoadData(const std::string& filepath);
	SimulationBuilder& SetIPC(int ipc);

	//Pattern methods
	Simulation Build();

};

