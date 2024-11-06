#pragma once

#include <vector>
#include <array>

struct ElementUniv
{
	std::vector<std::array<double, 4>> CsiDerivative{};
	std::vector<std::array<double, 4>> EtaDerivative{};
};