#include "ElementUniv.h"
#include "../utils/Gauss.h"

void ElementUniv::calcPcN()
{
	m_PcN = Matrix(m_IPC, 4);
	auto integrationPoints = Gauss::GetIntegrationPoints2D(m_IPC);

	for (size_t i = 0; i < m_IPC; ++i)
	{
		double csi = integrationPoints[i].Node.csi;
		double eta = integrationPoints[i].Node.eta;

		m_PcN(i, 0) = 0.25 * (1 - csi) * (1 - eta);
		m_PcN(i, 1) = 0.25 * (1 + csi) * (1 - eta);
		m_PcN(i, 2) = 0.25 * (1 + csi) * (1 + eta);
		m_PcN(i, 3) = 0.25 * (1 - csi) * (1 + eta);
	}
}

void ElementUniv::calcDerivatives()
{
	m_CsiDerivatives = Matrix(m_IPC, 4);
	m_EtaDerivatives = Matrix(m_IPC, 4);
	auto& integrationPoints = Gauss::GetIntegrationPoints2D(m_IPC);

	for (size_t i = 0; i < m_IPC; ++i)
	{
		double csi = integrationPoints[i].Node.csi;
		double eta = integrationPoints[i].Node.eta;

		m_CsiDerivatives(i, 0) = -0.25 * (1 - eta);
		m_CsiDerivatives(i, 1) = 0.25 * (1 - eta);
		m_CsiDerivatives(i, 2) = 0.25 * (1 + eta);
		m_CsiDerivatives(i, 3) = -0.25 * (1 + eta);

		m_EtaDerivatives(i, 0) = -0.25 * (1 - csi);
		m_EtaDerivatives(i, 1) = -0.25 * (1 + csi);
		m_EtaDerivatives(i, 2) = 0.25 * (1 + csi);
		m_EtaDerivatives(i, 3) = 0.25 * (1 - csi);
	}
}

ElementUniv::ElementUniv()
{
}

void ElementUniv::Initalize(const int ipc)
{
	if (m_Initialized)
		throw std::runtime_error("ERROR!: ElementUniv has already been initialized!");

	m_IPC = ipc;

	this->calcDerivatives();
	this->calcPcN();

	m_Initialized = true;
}

ElementUniv& ElementUniv::getInstance()
{
	static ElementUniv instance;
	return instance;
}

const Matrix& ElementUniv::getCsiDerivatives() const
{
	return m_CsiDerivatives;
}

const Matrix& ElementUniv::getEtaDerivatives() const
{
	return m_EtaDerivatives;
}

const Matrix& ElementUniv::getPcN() const
{
	return m_PcN;
}
