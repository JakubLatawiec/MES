#include "Matrix.h"

#include <iostream>
#include <stdexcept>
#include <iomanip>

Matrix Matrix::getMinor(size_t excludeRow, size_t excludeCol) const
{
	Matrix minor(m_Rows - 1, m_Cols - 1);
	size_t mRow = 0;
	size_t mCol = 0;
	for (size_t i = 0; i < m_Rows; ++i)
	{
		if (i == excludeRow)
			continue;
		mCol = 0;
		for (size_t j = 0; j < m_Cols; ++j)
		{
			if (j == excludeCol)
				continue;
			minor(mRow, mCol) = m_Data[i][j];
			++mCol;
		}
		++mRow;
	}
	return minor;
}

double& Matrix::operator()(size_t row, size_t col)
{
	if (row >= m_Rows)
		throw std::out_of_range("Row index out of range!");

	if (col >= m_Cols)
		throw std::out_of_range("Column index out of range!");

	return m_Data[row][col];
}

double Matrix::operator()(size_t row, size_t col) const
{
	if (row >= m_Rows)
		throw std::out_of_range("Row index out of range!");

	if (col >= m_Cols)
		throw std::out_of_range("Column index out of range!");

	return m_Data[row][col];
}

Matrix Matrix::operator+(const Matrix& other) const
{
	if (m_Rows != other.m_Rows || m_Cols != other.m_Cols)
		throw std::invalid_argument("MATRIX ADDITION: Incompatible matrix dimensions!");

	Matrix result(m_Rows, m_Cols);
	for (size_t i = 0; i < m_Rows; ++i)
		for (size_t j = 0; j < m_Cols; ++j)
			result(i, j) = m_Data[i][j] + other(i, j);

	return result;
}

Matrix Matrix::operator-(const Matrix& other) const
{
	if (m_Rows != other.m_Rows || m_Cols != other.m_Cols)
		throw std::invalid_argument("MATRIX SUBSTRACTION: Incompatible matrix dimensions!");

	Matrix result(m_Rows, m_Cols);
	for (size_t i = 0; i < m_Rows; ++i)
		for (size_t j = 0; j < m_Cols; ++j)
			result(i, j) = m_Data[i][j] - other(i, j);

	return result;
}

Matrix Matrix::operator*(const Matrix& other) const
{
	if (m_Cols != other.m_Rows)
		throw std::invalid_argument("MATRIX MULTIPLICATION: Incompatible matrix dimensions!");

	Matrix result(m_Rows, other.m_Rows);
	for (size_t i = 0; i < m_Rows; ++i)
		for (size_t j = 0; j < other.m_Cols; ++j)
			for (size_t k = 0; k < m_Cols; ++k)
				result(i, j) += m_Data[i][k] * other(k, j);

	return result;
}

Matrix Matrix::Transpose() const
{
	Matrix result(m_Cols, m_Rows);
	for (size_t i = 0; i < m_Rows; ++i)
		for (size_t j = 0; j < m_Cols; ++j)
			result(j, i) = m_Data[i][j];

	return result;
}

double Matrix::Determinant() const
{
	if (m_Rows != m_Cols)
		throw std::invalid_argument("MATRIX DETERMINANT: The matrix is ​​not a square matrix!");

	//Determinant for 1x1 matrix
	if (m_Rows == 1)
		return m_Data[0][0];
	//Determinant for 2x2 matrix
	else if (m_Rows == 2)
		return m_Data[0][0] * m_Data[1][1] - m_Data[0][1] * m_Data[1][0];
	//Deeterminant for NxN matrix (N > 2)
	else
	{
		double det = 0.0;
		for (size_t j = 0; j < m_Cols; ++j)
		{
			Matrix minor = getMinor(0, j);
			det += ((j % 2 == 0) ? 1 : -1) * m_Data[0][j] * minor.Determinant();
		}
		return det;
	}
}

Matrix Matrix::Inverse() const
{
	double det = Determinant();
	
	if (det == 0)
		throw std::runtime_error("MATRIX INVERSE: Determinant of matrix equals 0!");

	auto cofactor = [&]() -> Matrix
	{
		Matrix cof(m_Rows, m_Cols);
		for (size_t i = 0; i < m_Rows; ++i)
			for (size_t j = 0; j < m_Cols; ++j)
			{
				Matrix minor = getMinor(i, j);
				cof(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * minor.Determinant();
			}

		return cof;
	};

	Matrix cof = cofactor();
	Matrix adj = cof.Transpose();

	Matrix inv(m_Rows, m_Cols);
	for (size_t i = 0; i < m_Rows; ++i)
		for (size_t j = 0; j < m_Cols; ++j)
			inv(i, j) = adj(i, j) / det;

	return inv;
}

void Matrix::Display(int precision)
{
	std::cout << std::setprecision(precision);
	for (const auto& row : m_Data)
	{
		for (const auto& val : row)
			std::cout << val << "\t";
		std::cout << "\n";
	}
		
}

void Matrix::setCols(size_t size)
{
	for (auto& row : m_Data)
		row.resize(size, 0.0);
	m_Cols = size;
}

void Matrix::setRows(size_t size)
{
	m_Data.resize(size, std::vector<double>(m_Cols, 0.0));
	m_Rows = size;
}
