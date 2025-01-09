#pragma once

#include "../utils/Matrix.h"
#include "../utils/Gauss.h"
#include <array>
#include <iostream>


class ElementUniv
{
private:
	int m_IPC{};
	int m_IPC2D{};
	int m_ENC{};
	Matrix m_CsiDerivatives{};
	Matrix m_EtaDerivatives{};
	Matrix m_ShapeFunctionValues{};
	std::array<Matrix, 4> m_PVectors{};
	std::array<Matrix, 4> m_HbcMatrixes{};

	void calcDerivatives();
	void calcShapeFunctionValues();
	void calcSurface();

public:
	ElementUniv(int ipc, int enc = 4);
	
	inline const Matrix& getCsiDerivatives() const { return m_CsiDerivatives; }
	inline const Matrix& getEtaDerivatives() const { return m_EtaDerivatives; }
	inline const Matrix& getShapeFunctionValues() const { return m_ShapeFunctionValues; }
	inline const std::array<Matrix, 4>& getPVectors() const { return m_PVectors; }
	inline const std::array<Matrix, 4>& getHbcMatrixes() const { return m_HbcMatrixes; }
};