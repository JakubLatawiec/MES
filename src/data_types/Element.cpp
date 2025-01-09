#include "Element.h"

#include <iostream>
#include "CalculationsData.h"

int Element::m_IPC;
int Element::m_SFC;
std::vector<Coefficient2D> Element::m_IntegrationPoints;
Surface Element::m_Surface;

void Element::CalcSurface(int ipc)
{
	m_Surface.calcSurface(ipc);
}

void Element::CalcJacobians(const std::vector<Node>& nodes)
{
	//To change:
	//m_IntegrationPoints = Gauss::GetIntegrationPoints2D(m_IPC);

	const auto& ipc = CalculationsData::getInstance().getIPC2D();
	const auto& integrationPoints = CalculationsData::getInstance().getIntegrationPoints2D();

	const auto& csiDerivatives = ElementUniv::getInstance().getCsiDerivatives();
	const auto& etaDerivatives = ElementUniv::getInstance().getEtaDerivatives();
	m_JacobianMatrixes.resize(ipc);

	for (int i = 0; i < ipc; ++i)
	{
		Jacobian jacobian;
		for (size_t j = 0; j < 4; ++j)
		{
			int node_id = m_NodesID[j] - 1;
			jacobian.J(0, 0) += csiDerivatives(i, j) * nodes[j].x;
			jacobian.J(0, 1) += csiDerivatives(i, j) * nodes[j].y;
			jacobian.J(1, 0) += etaDerivatives(i, j) * nodes[j].x;
			jacobian.J(1, 1) += etaDerivatives(i, j) * nodes[j].y;
		}

		jacobian.DetJ = jacobian.J.Determinant();
		jacobian.J1 = jacobian.J.Inverse();

		m_JacobianMatrixes[i] = jacobian;
	}
}

void Element::CalcStiffnessMatrixes(double conductivity)
{
	const auto& ipc = CalculationsData::getInstance().getIPC2D();
	const auto& sfc = CalculationsData::getInstance().getSFC();
	const auto& integrationPoints = CalculationsData::getInstance().getIntegrationPoints2D();

	const auto& csiDerivatives = ElementUniv::getInstance().getCsiDerivatives();
	const auto& etaDerivatives = ElementUniv::getInstance().getEtaDerivatives();

	m_StiffnessMatrixes.resize(ipc);

	Matrix dNdX(ipc, sfc);
	Matrix dNdY(ipc, sfc);

	for (size_t i = 0; i < ipc; ++i)
	{
		for (size_t j = 0; j < sfc; ++j)
		{
			dNdX(i, j) = m_JacobianMatrixes[i].J1(0, 0) * csiDerivatives(i, j) + m_JacobianMatrixes[i].J1(0, 1) * etaDerivatives(i, j);
			dNdY(i, j) = m_JacobianMatrixes[i].J1(1, 0) * csiDerivatives(i, j) + m_JacobianMatrixes[i].J1(1, 1) * etaDerivatives(i, j);
		}
	}

	for (size_t i = 0; i < ipc; ++i)
		m_StiffnessMatrixes[i] = m_JacobianMatrixes[i].DetJ * conductivity * (dNdX.getRow(i).Transpose() * dNdX.getRow(i) + dNdY.getRow(i).Transpose() * dNdY.getRow(i));


		
	m_StiffnessMatrix = Matrix(sfc, sfc);
	for (size_t i = 0; i < ipc; ++i)
		m_StiffnessMatrix = m_StiffnessMatrix + m_StiffnessMatrixes[i] * m_IntegrationPoints[i].SurfArea;
}

void Element::CalcHbcMatrixes(const std::vector<Node>& nodes, double alpha)
{
	for (size_t side = 0; side < nodes.size(); ++side)
	{
		const Node& current = nodes[side];
		const Node& next = nodes[(side + 1) % nodes.size()];

		if (current.isBorderCondition && next.isBorderCondition)
		{
			double sideLength = std::sqrt(std::pow(next.x - current.x, 2) + std::powf(next.y - current.y, 2));
			double detJ = sideLength / 2.0;
			m_HbcMatrixes[side] = m_Surface.HbcMatrixes[side] * alpha * detJ;
		}
		else
			m_HbcMatrixes[side] = Matrix(4, 4);
	}

	m_StiffnessMatrix.Display();
	for (size_t side = 0; side < 4; ++side)
	{
		m_StiffnessMatrix = m_StiffnessMatrix + m_HbcMatrixes[side];
	}
}

void Element::CalcPVector(const std::vector<Node>& nodes, double alpha, double Tot)
{
	m_PVector = Matrix(4, 1);

	for (size_t side = 0; side < nodes.size(); ++side)
	{
		const Node& current = nodes[side];
		const Node& next = nodes[(side + 1) % nodes.size()];

		if (current.isBorderCondition && next.isBorderCondition)
		{
			double sideLength = std::sqrt(std::pow(next.x - current.x, 2) + std::powf(next.y - current.y, 2));
			double detJ = sideLength / 2.0;
			m_PVector = m_PVector + (m_Surface.PVectors[side] * alpha * Tot * detJ);
		}
	}
}

void Element::CalcCMatrix(double c, double rho)
{
	auto integrationPoints = Gauss::GetIntegrationPoints2D(m_IPC);
	m_CMatrixes.resize(m_IPC);

	auto& PcN = ElementUniv::getInstance().getPcN();

	for (size_t i = 0; i < m_IPC; ++i)
	{
		m_CMatrixes[i] = c * rho * (PcN.getRow(i).Transpose() * PcN.getRow(i)) * m_JacobianMatrixes[i].DetJ;
	}

	m_CMatrix = Matrix(4, 4);
	for (size_t i = 0; i < m_IPC; ++i)
	{
		m_CMatrix = m_CMatrix + m_CMatrixes[i] * integrationPoints[i].SurfArea;
	}
}


void Element::PrintJacobianMatrixes()
{
	std::cout << "JACOBIAN MATRIXES\n";
	for (size_t i = 0; i < m_IPC; ++i)
	{
		std::cout << "pc" << i + 1 << ":\n";
		m_JacobianMatrixes[i].J1.Display(6);
		std::cout << m_JacobianMatrixes[i].DetJ << "\n";
	}
}

void Element::PrintStiffnessMatrixes()
{
	std::cout << "STIFFNESS MATRIXES\n";
	for (size_t i = 0; i < m_IPC; ++i)
	{
		std::cout << "pc" << i + 1 << ":\n";
		m_StiffnessMatrixes[i].Display(6);
	}
}

void Element::PrintStiffnessMatrix()
{
	std::cout << "ELEMENT STIFFNESS MATRIX\n";
	m_StiffnessMatrix.Display(6);
}

void Element::PrintNodesID()
{
	std::cout << "NODES IDS\n";
	for (auto& nodeID : this->m_NodesID)
	{
		std::cout << nodeID << ", ";
	}
	std::cout << std::endl;
}

void Element::PrintHbcMatrixes()
{
	std::cout << "HBC MATRIXES\n";
	for (size_t i = 0; i < m_HbcMatrixes.size(); ++i)
	{
		m_HbcMatrixes[i].Display();
	}
}

void Element::PrintPVector()
{
	std::cout << "P VECTOR\n";
	m_PVector.Display();
}

void Element::PrintCMatrixes()
{
	std::cout << "C MATRIXES\n";
	for (size_t i = 0; i < m_CMatrixes.size(); ++i)
	{
		m_CMatrixes[i].Display();
		std::cout << "\n";
	}
}

void Element::PrintCMatrix()
{
	std::cout << "C MATRIX\n";
	m_CMatrix.Display();
}

void Element::setNodesID(size_t index, int id)
{
	m_NodesID[index] = id;
}

const std::array<int, 4>& Element::getNodesID() const
{
	return m_NodesID;
}

const Matrix& Element::getStifnessMatrix() const
{
	return m_StiffnessMatrix;
}

const Matrix& Element::getPVector() const
{
	return m_PVector;
}

const Matrix& Element::getCMatrix() const
{
	return m_CMatrix;
}
