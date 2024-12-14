#include "simulation_builder/SimulationBuilder.h"

#include <iostream>

int main()
{
	try
	{
		Simulation simulation = SimulationBuilder()
			.LoadData("../data/testDataMixGrid.txt")
			.SetIPC(9)
			.SetSurfaceIPC(2)
			.Build();

		simulation.Run();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << "\n";
	}

	return 0;
}
