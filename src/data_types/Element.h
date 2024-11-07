#pragma once

#include <array>
#include <vector>

#include "Jacobian.h"
#include "ElementUniv.h"
#include "Node.h"
#include "../utils/Gauss.h"

class Element
{
private:
	//Static variables
	static ElementUniv m_ElementUniv;
	static int m_IPC;
	static int m_SFC;
	static std::vector<Coefficient2D> m_IntegrationPoints;

private:
	//Object variables
	std::array<int, 4> m_NodesID{};
	std::vector<Jacobian> m_JacobianMatrixes{};
	std::vector<Matrix> m_StiffnessMatrixes{};
	Matrix m_StiffnessMatrix{};


public:
	//Constructors
	Element() = default;
	
	//Calculation methods
	static void CalcElementUniv(int ipc);
	void CalcJacobians(const std::vector<Node>& nodes);
	void CalcStiffnessMatrixes(double conductivity);

	//Debug methods
	static void PrintElementUniv();
	void PrintJacobianMatrixes();
	void PrintStiffnessMatrixes();
	void PrintStiffnessMatrix();

	//Setters
	void setNodesID(size_t index, int id);

	//Getters
	const std::array<int, 4>& getNodesID() const;
};