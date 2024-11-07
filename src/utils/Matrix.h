#pragma once

#include <vector>

class Matrix
{
private:
	std::vector<std::vector<double>> m_Data{};
	size_t m_Rows{};
	size_t m_Cols{};

	Matrix getMinor(size_t excludeRow, size_t excludeCol) const;

public:
	Matrix(size_t rows, size_t cols)
		: m_Rows(rows), m_Cols(cols), m_Data(rows, std::vector<double>(cols, 0.0)) {}

	Matrix() = default;

	//Read and write operator
	double& operator()(size_t row, size_t col);

	//Read only operator
	double operator()(size_t row, size_t col) const;

	//Addition operator
	Matrix operator+(const Matrix& other) const;

	//Substraction operator
	Matrix operator-(const Matrix& other) const;

	//Multiplication operator
	Matrix operator*(const Matrix& other) const;

	//Transpose method
	Matrix Transpose() const;

	//Calculating determinant
	double Determinant() const;

	//Inversing method
	Matrix Inverse() const;

	//Debug methods
	void Display(int precision = 8);

	//Setters
	void setCols(size_t size);
	void setRows(size_t size);

	//Getters
	inline size_t getCols() const { return m_Cols; }
	inline size_t getRows() const { return m_Rows; }
};

