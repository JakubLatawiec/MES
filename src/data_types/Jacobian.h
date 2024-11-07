#pragma once

#include "../utils/Matrix.h"

struct Jacobian
{
	Matrix J;
	Matrix J1;
	double DetJ{};

	Jacobian()
		: J(Matrix(2,2)), J1(Matrix(2,2)) {}
};