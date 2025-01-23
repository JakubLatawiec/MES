#include "SimulationBuilder.h"

#include <stdexcept>

#include "../data_parser/DataParser.h"

SimulationBuilder& SimulationBuilder::LoadData(const std::string& filepath)
{
	DataParser::ParseData(filepath, m_Simulation.m_GlobalData, m_Simulation.m_Grid);
	return *this;
}

SimulationBuilder& SimulationBuilder::SetIPC(int ipc)
{
	double sqrtValue = sqrt(ipc);
	if (floor(sqrtValue) != sqrtValue)
		throw std::invalid_argument("SIMULATION BUILDER: Invalid integration points count!");

	m_Simulation.m_IPC = ipc;
	return *this;
}

SimulationBuilder& SimulationBuilder::SetSurfaceIPC(int ipc)
{
	m_Simulation.m_SurfaceIPC = ipc;
	return *this;
}

SimulationBuilder& SimulationBuilder::UseParaView(bool useParaView)
{
	m_Simulation.m_UseParaView = useParaView;
	return *this;
}

Simulation SimulationBuilder::Build()
{
	return m_Simulation;
}
