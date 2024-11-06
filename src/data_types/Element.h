#pragma once

#include <array>
#include <vector>

#include "Jacobian.h"

struct Element
{
	std::array<int, 4> NodesID{};
	std::vector<Jacobian> Jacobians{};
};