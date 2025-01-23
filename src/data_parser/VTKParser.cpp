#include "VTKParser.h"

#include <fstream>
#include <iomanip>

void VTKParser::parseToVTK(const Grid& grid)
{
	std::ofstream outputFile("../data/output/model.vtk");
	if (!outputFile.is_open()) {
		std::cerr << "Error: Cannot create output file." << std::endl;
		return;
	}

	outputFile << "# vtk DataFile Version 3.0\n";
	outputFile << "Finite Element Mesh\n";
	outputFile << "ASCII\n";
	outputFile << "DATASET UNSTRUCTURED_GRID\n";

	const auto& nodes = grid.Nodes;
	outputFile << "POINTS " << nodes.size() << " float\n";
	for (const auto& node : nodes) {
		outputFile << std::fixed << std::setprecision(6) << node.x << " " << node.y << " 0.0\n";
	}

    const auto& elements = grid.Elements;
    size_t totalIndices = 0;
    for (const auto& element : elements) {
        totalIndices += element.getNodesID().size() + 1;
    }
    outputFile << "CELLS " << elements.size() << " " << totalIndices << "\n";
    for (const auto& element : elements) {
        const auto& nodeIDs = element.getNodesID();
        outputFile << nodeIDs.size();
        for (const auto& nodeID : nodeIDs) {
            outputFile << " " << (nodeID - 1);
        }
        outputFile << "\n";
    }

    outputFile << "CELL_TYPES " << elements.size() << "\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        outputFile << "9\n";
    }

    outputFile.close();
}

void VTKParser::SaveTemperaturesToVTK(const std::string& filename, const Grid& grid, const Matrix& temperatures)
{
    std::ofstream outputFile("../data/output/" + filename + ".vtk");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Cannot create output file." << std::endl;
        return;
    }


    outputFile << "# vtk DataFile Version 3.0\n";
    outputFile << "Finite Element Mesh\n";
    outputFile << "ASCII\n";
    outputFile << "DATASET UNSTRUCTURED_GRID\n";

    const auto& nodes = grid.Nodes;
    outputFile << "POINTS " << nodes.size() << " float\n";
    for (const auto& node : nodes) {
        outputFile << std::fixed << std::setprecision(6) << node.x << " " << node.y << " 0.0\n";
    }

    const auto& elements = grid.Elements;
    size_t totalIndices = 0;
    for (const auto& element : elements) {
        totalIndices += element.getNodesID().size() + 1;
    }
    outputFile << "CELLS " << elements.size() << " " << totalIndices << "\n";
    for (const auto& element : elements) {
        const auto& nodeIDs = element.getNodesID();
        outputFile << nodeIDs.size();
        for (const auto& nodeID : nodeIDs) {
            outputFile << " " << (nodeID - 1);
        }
        outputFile << "\n";
    }

    outputFile << "CELL_TYPES " << elements.size() << "\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        outputFile << "9\n";
    }

    outputFile << "POINT_DATA " << temperatures.getRowsSize() << "\n";
    outputFile << "SCALARS temperature float 1\n";
    outputFile << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < temperatures.getRowsSize(); ++i) {
        outputFile << std::fixed << std::setprecision(6) << temperatures(i, 0) << "\n";
    }

    outputFile.close();
}

void VTKParser::WritePVDFile(size_t numSteps, double timeStep)
{
    std::ofstream outputFile("../data/output/simulation.pvd");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Cannot create output file." << std::endl;
        return;
    }

    outputFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    outputFile << "  <Collection>\n";

    // Write datasets for each time step
    for (size_t step = 0; step < numSteps; ++step) {
        double time = step * timeStep;
        outputFile << "    <DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"temperatures_step_" << step << ".vtk\"/>\n";
    }

    outputFile << "  </Collection>\n";
    outputFile << "</VTKFile>\n";

    outputFile.close();
}
