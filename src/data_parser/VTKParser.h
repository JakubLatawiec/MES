#pragma once

#include "../data_types/Grid.h"

class VTKParser
{
public:
	static void parseToVTK(const Grid& grid);
	static void SaveTemperaturesToVTK(const std::string& filename, const Grid& grid, const Matrix& temperatures);
	static void WritePVDFile(size_t numSteps, double timeStep);
};

