#include "Element.h"

#include <iostream>

#include "../utils/Gauss.h"

ElementUniv Element::m_ElementUniv;

void Element::CalcJacobians(const std::vector<Node>& nodes)
{
	std::cout << "CALC JACOBIANS:\n";

	int npc = m_ElementUniv.EtaDerivative.size();

	this->Jacobians.resize(npc);

	for (int i = 0; i < npc; ++i)
	{
		Jacobian jacobian;
		for (int j = 0; j < 4; ++j)
		{
			int node_id = this->NodesID[j] - 1;
			jacobian.J[0][0] += m_ElementUniv.CsiDerivative[i][j] * nodes[node_id].x;
			jacobian.J[0][1] += m_ElementUniv.CsiDerivative[i][j] * nodes[node_id].y;
			jacobian.J[1][0] += m_ElementUniv.EtaDerivative[i][j] * nodes[node_id].x;
			jacobian.J[1][1] += m_ElementUniv.EtaDerivative[i][j] * nodes[node_id].y;
		}

		jacobian.DetJ = jacobian.J[0][0] * jacobian.J[1][1] - jacobian.J[0][1] * jacobian.J[1][0];

		jacobian.J1[0][0] = jacobian.J[1][1] * (1.0 / jacobian.DetJ);
		jacobian.J1[0][1] = -jacobian.J[0][1] * (1.0 / jacobian.DetJ);
		jacobian.J1[1][0] = -jacobian.J[1][0] * (1.0 / jacobian.DetJ);
		jacobian.J1[1][1] = jacobian.J[0][0] * (1.0 / jacobian.DetJ);

		this->Jacobians[i] = jacobian;

		std::cout << this->Jacobians[i].J[0][0] << " " << this->Jacobians[i].J[0][1] << "\n";
		std::cout << this->Jacobians[i].J[1][0] << " " << this->Jacobians[i].J[1][1] << "\n";
		std::cout << this->Jacobians[i].J1[0][0] << " " << this->Jacobians[i].J1[0][1] << "\n";
		std::cout << this->Jacobians[i].J1[1][0] << " " << this->Jacobians[i].J1[1][1] << "\n";
	}

	
}

void Element::CalcElementUniv(int npc)
{
	m_ElementUniv.CsiDerivative.resize(npc);
	m_ElementUniv.EtaDerivative.resize(npc);

	auto integrationPoints = Gauss::GetIntegrationPoints2D(npc);

	for (int i = 0; i < npc; ++i)
	{
		m_ElementUniv.CsiDerivative[i][0] = -0.25 * (1 - integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative[i][1] =  0.25 * (1 - integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative[i][2] =  0.25 * (1 + integrationPoints[i].Node.eta);
		m_ElementUniv.CsiDerivative[i][3] = -0.25 * (1 + integrationPoints[i].Node.eta);

		m_ElementUniv.EtaDerivative[i][0] = -0.25 * (1 - integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative[i][1] = -0.25 * (1 + integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative[i][2] =  0.25 * (1 + integrationPoints[i].Node.csi);
		m_ElementUniv.EtaDerivative[i][3] =  0.25 * (1 - integrationPoints[i].Node.csi);
	}
}

void Element::PrintElementUniv()
{
	std::cout << "ELEMENT UNIV: CsiDerivative\n";
	for (auto& pc : m_ElementUniv.CsiDerivative)
	{
		for (auto& e : pc)
			std::cout << e << "\t";
		std::cout << std::endl;
	}

	std::cout << "ELEMENT UNIV: EtaDerivative\n";
	for (auto& pc : m_ElementUniv.EtaDerivative)
	{
		for (auto& e : pc)
			std::cout << e << "\t";
		std::cout << std::endl;
	}
}
