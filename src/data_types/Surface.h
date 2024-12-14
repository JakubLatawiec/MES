#pragma once

#include <array>
#include "../utils/Matrix.h"

struct Surface
{
	std::array<Matrix, 4> PVectors;
	std::array<Matrix, 4> HbcMatrixes;
	void calcSurface(int npc);
};