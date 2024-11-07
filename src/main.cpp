#include "simulation_builder/SimulationBuilder.h"

#include <iostream>

int main()
{
	try
	{
		Simulation simulation = SimulationBuilder()
			.LoadData("../data/testDataSimplex2.txt")
			.SetIPC(4)
			.Build();

		simulation.Run();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << "\n";
	}

	return 0;
}
