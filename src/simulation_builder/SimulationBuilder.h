#pragma once

#include "../data_parser/DataParser.h"
#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

#include "Simulation.h"

class SimulationBuilder
{
private:
	Simulation simulation;

public:
	SimulationBuilder& LoadData(const std::string& filepath);

	Simulation Build();

};

