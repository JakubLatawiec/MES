#pragma once

#include <vector>
#include <array>

#include "../utils/Matrix.h"

struct ElementUniv
{
	//std::vector<std::array<double, 4>> CsiDerivative{};
	//std::vector<std::array<double, 4>> EtaDerivative{};
	Matrix CsiDerivative;
	Matrix EtaDerivative;
};