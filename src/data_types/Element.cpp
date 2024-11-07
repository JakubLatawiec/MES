#include "Element.h"

#include <iostream>

#include "../utils/Gauss.h"

ElementUniv Element::m_ElementUniv;

void Element::CalcJacobians(const std::vector<Node>& nodes)
{
	std::cout << "CALC JACOBIANS:\n";

	int npc = m_ElementUniv.CsiDerivative.getRows();

	//int npc = m_ElementUniv.EtaDerivative.size();

	this->Jacobians.resize(npc);

	for (int i = 0; i < npc; ++i)
	{
		Jacobian jacobian;
		for (int j = 0; j < nodes.size(); ++j)
		{
			int node_id = this->NodesID[j] - 1;
			jacobian.J(0, 0) += m_ElementUniv.CsiDerivative(i, j) * nodes[node_id].x;
			jacobian.J(0, 1) += m_ElementUniv.CsiDerivative(i, j) * nodes[node_id].y;
			jacobian.J(1, 0) += m_ElementUniv.EtaDerivative(i, j) * nodes[node_id].x;
			jacobian.J(1, 1) += m_ElementUniv.EtaDerivative(i, j) * nodes[node_id].y;
		}

		jacobian.DetJ = jacobian.J.Determinant();
		jacobian.J1 = jacobian.J.Inverse();

		this->Jacobians[i] = jacobian;

		this->Jacobians[i].J1.Display();
		this->Jacobians[i].J.Display();
		std::cout << this->Jacobians[i].DetJ << "\n";
	}

	
}

void Element::CalcElementUniv(int npc)
{
	m_ElementUniv.CsiDerivative = Matrix(npc, 4);
	m_ElementUniv.EtaDerivative = Matrix(npc, 4);

	auto integrationPoints = Gauss::GetIntegrationPoints2D(npc);

	for (int i = 0; i < npc; ++i)
	{
		m_ElementUniv.CsiDerivative(i, 0) = -0.25 * (1 - integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 1) =  0.25 * (1 - integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 2) =  0.25 * (1 + integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative(i, 3) = -0.25 * (1 + integrationPoints[i].Node.eta);

		m_ElementUniv.EtaDerivative(i, 0) = -0.25 * (1 - integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 1) = -0.25 * (1 + integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 2) =  0.25 * (1 + integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative(i, 3) =  0.25 * (1 - integrationPoints[i].Node.csi);
	}
}

void Element::PrintElementUniv()
{
	std::cout << "ELEMENT UNIV: CsiDerivative\n";
	m_ElementUniv.CsiDerivative.Display(6);

	std::cout << "ELEMENT UNIV: EtaDerivative\n";
	m_ElementUniv.EtaDerivative.Display(6);
}
