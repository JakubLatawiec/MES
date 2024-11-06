#pragma once

#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

class Simulation
{
private:
	Simulation() = default;
	friend class SimulationBuilder;

	void printLoadedData();
	void calcElementJacobians();

public:
	GlobalData GlobalData{};
	Grid Grid{};
	void Run();
};

