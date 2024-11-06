#include <iostream>

#include "simulation_builder/SimulationBuilder.h"

int main()
{
	Simulation simulation = SimulationBuilder()
		.loadData("../data/testData.txt")
		.Build();

	simulation.Run();

	return 0;
}
