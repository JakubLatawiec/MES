#pragma once

#include "../utils/Matrix.h"
#include "../utils/Gauss.h"
#include <array>
#include <iostream>


class ElementUniv
{
private:
	int m_IPC{};
	bool m_Initialized = false;

	Matrix m_CsiDerivatives{};
	Matrix m_EtaDerivatives{};
	Matrix m_PcN{};
	std::array<Matrix, 4> m_PVectors{};
	std::array<Matrix, 4> m_HbcMatrixes{};

	void calcPcN();
	void calcDerivatives();

	ElementUniv();
	~ElementUniv() {};
	ElementUniv(const ElementUniv&) = delete;
	ElementUniv& operator=(const ElementUniv&) = delete;

public:
	void Initalize(const int ipc);
	static ElementUniv& getInstance();
	const Matrix& getCsiDerivatives() const;
	const Matrix& getEtaDerivatives() const;
	const Matrix& getPcN() const;
};