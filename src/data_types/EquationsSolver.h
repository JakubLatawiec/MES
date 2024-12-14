#pragma once

#include "../utils/Matrix.h"

class EquationsSolver
{
private:
	Matrix m_GlobalStiffnessMatrix{};
	Matrix m_GlobalPVector{};
	Matrix m_GlobalCMatrix{};
	Matrix m_L{};
	Matrix m_U{};
	Matrix m_GlobalTVector{};

	//Calculations methods
	void luDecomposition();
	Matrix forwardSubstitution();
	Matrix backwardSubstitution(Matrix& y);

	//Calc new
	std::pair<Matrix, Matrix> LUDecomposition(const Matrix& A);
	Matrix ForwardSubstitution(const Matrix& L, const Matrix& b);
	Matrix BackwardSubstitution(const Matrix& U, const Matrix& y);


public:
	EquationsSolver(int nodeCount);
	EquationsSolver() = default;

	//Calculations methods
	void SolveEquation(double dt, double simulationTime, double initailTemp);

	//Setters
	void setGlobalStiffnessMatrix(size_t row, size_t col, double val);
	void setGlobalPVector(size_t index, double val);
	void setGlobalCMatrix(size_t row, size_t col, double val);

	//Debug
	void PrintGlobalStiffnessMatrix();
	void PrintGlobalPVector();
	void PrintGlobalTVector();
	void PrintGlobalCMatrix();
};

