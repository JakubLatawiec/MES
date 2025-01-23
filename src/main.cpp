#include "simulation_builder/SimulationBuilder.h"

#include <iostream>

int main()
{
	try
	{
		Simulation simulation = SimulationBuilder()
			.LoadData("../data/wall.txt")
			.SetIPC(16)
			.SetSurfaceIPC(4)
			.UseParaView()
			.Build();

		simulation.Run();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << "\n";
	}

	return 0;
}
