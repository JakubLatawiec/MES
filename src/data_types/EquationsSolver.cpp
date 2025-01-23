#include "EquationsSolver.h"
#include <stdexcept>
#include <iostream>

std::pair<Matrix, Matrix> EquationsSolver::LUDecomposition(const Matrix& A)
{
	size_t n = A.getRowsSize();
	if (n != A.getColsSize())
		throw std::invalid_argument("LU DECOMPOSITION: Matrix A is not square matrix!");

	Matrix L = Matrix(n, n);
	Matrix U = Matrix(A);

	for (size_t i = 0; i < n; ++i)
	{
		L(i, i) = 1.0;

		for (size_t j = i + 1; j < n; ++j)
		{
			if (U(i, i) == 0)
				throw std::runtime_error("LU DECOMPOSITION: Division by 0!");
			double factor = U(j, i) / U(i, i);
			L(j, i) = factor;

			for (size_t k = i; k < n; ++k)
				U(j, k) -= factor * U(i, k);
		}
	}

	return { L, U };
}

Matrix EquationsSolver::ForwardSubstitution(const Matrix& L, const Matrix& b)
{
	size_t n = L.getRowsSize();
	Matrix y = Matrix(n, 1);

	for (size_t i = 0; i < n; ++i)
	{
		double sum = 0.0;
		for (size_t j = 0; j < i; ++j)
			sum += L(i, j) * y(j, 0);

		if (L(i, i) == 0)
			throw std::runtime_error("FORWARD SUBSTITUTION: Division by 0!");
		y(i, 0) = (b(i, 0) - sum) / L(i, i);
	}

	return y;
}

Matrix EquationsSolver::BackwardSubstitution(const Matrix& U, const Matrix& y)
{
	size_t n = U.getRowsSize();
	Matrix x = Matrix(n, 1);

	for (int i = n - 1; i >= 0; --i)
	{
		double sum = 0.0;
		for (size_t j = i + 1; j < n; ++j)
			sum += U(i, j) * x(j, 0);

		if (U(i, i) == 0)
			throw std::runtime_error("BACKWARD SUBSTITUTION: Division by 0!");

		x(i, 0) = (y(i, 0) - sum) / U(i, i);
	}

	return x;
}

EquationsSolver::EquationsSolver(int nodeCount)
{
	m_GlobalStiffnessMatrix = Matrix(nodeCount, nodeCount);
	m_GlobalPVector = Matrix(nodeCount, 1);
	m_GlobalCMatrix = Matrix(nodeCount, nodeCount);
}

std::vector<std::pair<double, Matrix>> EquationsSolver::SolveEquation(double dt, double simulationTime, double initialTemp)
{
	Matrix A = m_GlobalStiffnessMatrix + (m_GlobalCMatrix * (1.0 / dt));
	auto [L, U] = LUDecomposition(A);

	Matrix t_prev = Matrix(A.getRowsSize(), 1);
	for (size_t i = 0; i < A.getRowsSize(); ++i)
		t_prev(i, 0) = initialTemp;

	std::vector<std::pair<double, Matrix>> results;

	results.emplace_back(0.0, t_prev);

	for (double currentTime = dt; currentTime <= simulationTime; currentTime += dt)
	{
		Matrix b = (m_GlobalCMatrix * (1.0 / dt)) * t_prev + m_GlobalPVector;
		Matrix y = ForwardSubstitution(L, b);
		Matrix x = BackwardSubstitution(U, y);

		results.emplace_back(currentTime, x);

		t_prev = x;
	}

	return results;
}

void EquationsSolver::setGlobalStiffnessMatrix(size_t row, size_t col, double val)
{
	m_GlobalStiffnessMatrix(row, col) += val;
}

void EquationsSolver::setGlobalPVector(size_t index, double val)
{
	m_GlobalPVector(index, 0) += val;
}

void EquationsSolver::setGlobalCMatrix(size_t row, size_t col, double val)
{
	m_GlobalCMatrix(row, col) += val;
}

void EquationsSolver::PrintGlobalStiffnessMatrix()
{
	m_GlobalStiffnessMatrix.Display();
}

void EquationsSolver::PrintGlobalPVector()
{
	m_GlobalPVector.Display();
}

void EquationsSolver::PrintGlobalTVector()
{
	m_GlobalTVector.Display();
}

void EquationsSolver::PrintGlobalCMatrix()
{
	m_GlobalCMatrix.Display();
}
