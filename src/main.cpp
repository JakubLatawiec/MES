#include <iostream>

#include "simulation_builder/SimulationBuilder.h"

int main()
{
	Simulation simulation = SimulationBuilder()
		.LoadData("../data/testDataSimplex.txt")
		.Build();

	simulation.Run();

	return 0;
}
