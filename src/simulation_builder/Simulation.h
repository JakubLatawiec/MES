#pragma once

#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

class Simulation
{
private:
	Simulation() = default;
	friend class SimulationBuilder;

	void printLoadedData();

public:
	GlobalData GlobalData{};
	Grid Grid{};
	void Run();
};

