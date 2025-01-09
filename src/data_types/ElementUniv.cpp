#include "ElementUniv.h"
#include "../utils/Gauss.h"

void ElementUniv::calcDerivatives()
{
	m_CsiDerivatives = Matrix(m_IPC2D, m_ENC);
	m_EtaDerivatives = Matrix(m_IPC2D, m_ENC);
	auto& integrationPoints = Gauss::getIntegrationPoints2D();

	for (size_t i = 0; i < m_IPC2D; ++i)
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

void ElementUniv::calcShapeFunctionValues()
{
	m_ShapeFunctionValues = Matrix(m_IPC2D, m_ENC);
	auto& integrationPoints = Gauss::getIntegrationPoints2D();

	for (size_t i = 0; i < m_IPC2D; ++i)
	{
		double csi = integrationPoints[i].Node.csi;
		double eta = integrationPoints[i].Node.eta;

		m_ShapeFunctionValues(i, 0) = 0.25 * (1 - csi) * (1 - eta);
		m_ShapeFunctionValues(i, 1) = 0.25 * (1 + csi) * (1 - eta);
		m_ShapeFunctionValues(i, 2) = 0.25 * (1 + csi) * (1 + eta);
		m_ShapeFunctionValues(i, 3) = 0.25 * (1 - csi) * (1 + eta);
	}
}

void ElementUniv::calcSurface()
{
	for (int side = 0; side < 4; ++side)
	{
		double csi = 0;
		double eta = 0;
		auto& integrationPoints = Gauss::getIntegrationPoints1D();

		m_HbcMatrixes[side] = Matrix(m_ENC, m_ENC);
		m_PVectors[side] = Matrix(m_ENC, 1);

		Matrix PcN = Matrix(m_IPC, m_ENC);

		for (int i = 0; i < m_IPC; ++i)
		{
			switch (side)
			{
			case 0:
				csi = integrationPoints[i].X;
				eta = -1;
				break;
			case 1:
				csi = 1;
				eta = integrationPoints[i].X;
				break;
			case 2:
				csi = integrationPoints[i].X;
				eta = 1;
				break;
			case 3:
				csi = -1;
				eta = integrationPoints[i].X;
				break;
			}

			PcN(i, 0) = 0.25 * (1 - csi) * (1 - eta);
			PcN(i, 1) = 0.25 * (1 + csi) * (1 - eta);
			PcN(i, 2) = 0.25 * (1 + csi) * (1 + eta);
			PcN(i, 3) = 0.25 * (1 - csi) * (1 + eta);
		}

		for (int i = 0; i < m_IPC; ++i)
		{
			m_HbcMatrixes[side] = m_HbcMatrixes[side] + (PcN.getRow(i).Transpose() * PcN.getRow(i) * integrationPoints[i].W);
			m_PVectors[side] = m_PVectors[side] + (PcN.getRow(i).Transpose() * integrationPoints[i].W);
		}
	}
}

ElementUniv::ElementUniv(int ipc, int enc) : m_IPC(ipc), m_ENC(enc)
{
	m_IPC2D = pow(ipc, 2);
	calcDerivatives();
	calcSurface();
	calcShapeFunctionValues();
}
