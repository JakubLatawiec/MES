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

public:
	std::array<int, 4> NodesID{};
	std::vector<Jacobian> Jacobians{};

	void CalcJacobians(const std::vector<Node>& nodes);

	Element() = default;

	static void CalcElementUniv(int npc);
	static void PrintElementUniv();

	//Getters
	inline static ElementUniv& GetElementUniv() { return m_ElementUniv; }
};