#include "SimulationBuilder.h"

#include "../data_parser/DataParser.h"

SimulationBuilder& SimulationBuilder::LoadData(const std::string& filepath)
{
	DataParser::ParseData(filepath, this->simulation.GlobalData, this->simulation.Grid);
	return *this;
}

Simulation SimulationBuilder::Build()
{
	return this->simulation;
}
