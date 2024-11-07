#pragma once

#include <array>
#include <vector>

#include "Jacobian.h"
#include "ElementUniv.h"
#include "Node.h"

class Element
{
private:
	static ElementUniv m_ElementUniv;
	std::vector<Matrix> m_StiffnessMatrixes;
	Matrix m_StiffnessMatrix;

public:
	std::array<int, 4> NodesID{};
	std::vector<Jacobian> Jacobians{};

	void CalcJacobians(const std::vector<Node>& nodes);
	void CalcStiffnessMatrixes(double conductivity);

	Element() = default;

	static void CalcElementUniv(int npc);
	static void PrintElementUniv();
};