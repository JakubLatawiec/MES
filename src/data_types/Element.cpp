#include "Element.h"

#include <iostream>

ElementUniv Element::m_ElementUniv;
int Element::m_IPC;
int Element::m_SFC;
std::vector<Coefficient2D> Element::m_IntegrationPoints;

void Element::CalcElementUniv(int ipc)
{
	m_IPC = ipc;
	m_SFC = 4;
	m_ElementUniv.CsiDerivative = Matrix(m_IPC, m_SFC);
	m_ElementUniv.EtaDerivative = Matrix(m_IPC, m_SFC);
	m_IntegrationPoints = Gauss::GetIntegrationPoints2D(ipc);

	for (size_t i = 0; i < ipc; ++i)
	{
		m_ElementUniv.CsiDerivative(i, 0) = -0.25 * (1 - m_IntegrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 1) = 0.25 * (1 - m_IntegrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 2) = 0.25 * (1 + m_IntegrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 3) = -0.25 * (1 + m_IntegrationPoints[i].Node.eta);

		m_ElementUniv.EtaDerivative(i, 0) = -0.25 * (1 - m_IntegrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 1) = -0.25 * (1 + m_IntegrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 2) = 0.25 * (1 + m_IntegrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 3) = 0.25 * (1 - m_IntegrationPoints[i].Node.csi);
	}
}

void Element::CalcJacobians(const std::vector<Node>& nodes)
{
	m_JacobianMatrixes.resize(m_IPC);
	for (int i = 0; i < m_IPC; ++i)
	{
		Jacobian jacobian;
		for (size_t j = 0; j < nodes.size(); ++j)
		{
			int node_id = m_NodesID[j] - 1;
			jacobian.J(0, 0) += m_ElementUniv.CsiDerivative(i, j) * nodes[j].x;
			jacobian.J(0, 1) += m_ElementUniv.CsiDerivative(i, j) * nodes[j].y;
			jacobian.J(1, 0) += m_ElementUniv.EtaDerivative(i, j) * nodes[j].x;
			jacobian.J(1, 1) += m_ElementUniv.EtaDerivative(i, j) * nodes[j].y;
		}

		jacobian.DetJ = jacobian.J.Determinant();
		jacobian.J1 = jacobian.J.Inverse();

		m_JacobianMatrixes[i] = jacobian;
	}
}

void Element::CalcStiffnessMatrixes(double conductivity)
{
	m_StiffnessMatrixes.resize(m_IPC);

	Matrix dNdX(m_IPC, m_SFC);
	Matrix dNdY(m_IPC, m_SFC);
	for (size_t i = 0; i < m_IPC; ++i)
	{
		for (size_t j = 0; j < m_SFC; ++j)
		{
			dNdX(i, j) = m_JacobianMatrixes[i].J1(0, 0) * m_ElementUniv.CsiDerivative(i, j) + m_JacobianMatrixes[i].J1(0, 1) * m_ElementUniv.EtaDerivative(i, j);
			dNdY(i, j) = m_JacobianMatrixes[i].J1(1, 0) * m_ElementUniv.CsiDerivative(i, j) + m_JacobianMatrixes[i].J1(1, 1) * m_ElementUniv.EtaDerivative(i, j);
		}
	}

	for (size_t i = 0; i < m_IPC; ++i)
		m_StiffnessMatrixes[i] = m_JacobianMatrixes[i].DetJ * conductivity * (dNdX.getRow(i).Transpose() * dNdX.getRow(i) + dNdY.getRow(i).Transpose() * dNdY.getRow(i));
		
	m_StiffnessMatrix = Matrix(m_SFC, m_SFC);
	for (size_t i = 0; i < m_IPC; ++i)
		m_StiffnessMatrix = m_StiffnessMatrix + m_StiffnessMatrixes[i] * m_IntegrationPoints[i].SurfArea;
}



void Element::PrintElementUniv()
{
	std::cout << "ELEMENT UNIV\n";
	std::cout << "CsiDerivative:\n";
	m_ElementUniv.CsiDerivative.Display(6);
	std::cout << "EtaDerivative\n";
	m_ElementUniv.EtaDerivative.Display(6);
	std::cout << std::endl;
}

void Element::PrintJacobianMatrixes()
{
	std::cout << "JACOBIAN MATRIXES\n";
	for (size_t i = 0; i < m_IPC; ++i)
	{
		std::cout << "pc" << i + 1 << ":\n";
		m_JacobianMatrixes[i].J1.Display(6);
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
	std::cout << "STIFFNESS MATRIXES\n";
	m_StiffnessMatrix.Display(6);
}

void Element::setNodesID(size_t index, int id)
{
	m_NodesID[index] = id;
}

const std::array<int, 4>& Element::getNodesID() const
{
	return m_NodesID;
}
