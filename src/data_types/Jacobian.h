#pragma once

#include <array>

struct Jacobian
{
	std::array<std::array<double, 2>, 2> J{};
	std::array<std::array<double, 2>, 2> J1{};
	double DetJ{};
};