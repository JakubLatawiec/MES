#include "Gauss.h"

#include <iostream>

std::unordered_map<int, std::vector<Coefficient1D>> Gauss::m_CoefficientsDatabase =
{
	{1, {
			{0, 2}
		}
	},
	{2, { 
			{-sqrt(1.0 / 3.0), 1.0}, 
			{sqrt(1.0 / 3.0), 1.0} 
		} 
	},
	{3, {
			{-sqrt(3.0 / 5.0), 5.0 / 9.0},
			{0.0, 8.0 / 9.0},
			{sqrt(3.0 / 5.0), 5.0 / 9.0}
		}
	},
	{4, {
			{-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), (18.0 - sqrt(30.0)) / 36.0},
			{-sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), (18.0 + sqrt(30.0)) / 36.0},
			{ sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), (18.0 + sqrt(30.0)) / 36.0},
			{ sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), (18.0 - sqrt(30.0)) / 36.0}
		}
	},
	{5, {
			{-(1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)), (322.0 - 13.0 * sqrt(70)) / 900.0},
			{-(1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)), (322.0 + 13.0 * sqrt(70)) / 900.0},
			{ 0.0, 128.0 / 225.0},
			{ (1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)), (322.0 + 13.0 * sqrt(70)) / 900.0},
			{ (1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)), (322.0 - 13.0 * sqrt(70)) / 900.0}
		}
	}
};

int Gauss::m_IPC{};
std::vector<Coefficient1D> Gauss::m_Coefficients1D{};
std::vector<Coefficient2D> Gauss::m_Coefficients2D{};

void Gauss::calcCoefficients1D()
{
	m_Coefficients1D.resize(m_IPC);
	m_Coefficients1D = m_CoefficientsDatabase[m_IPC];
}

void Gauss::calcCoefficients2D()
{
	m_Coefficients2D.resize(pow(m_IPC, 2));

	const auto& points = m_CoefficientsDatabase[m_IPC];
	for (const auto& csiPoint : points)
		for (const auto& etaPoint : points)
		{
			Coefficient2D coeff;
			coeff.Node.csi = csiPoint.X;
			coeff.Node.eta = etaPoint.X;
			coeff.SurfArea = csiPoint.W * etaPoint.W;
			m_Coefficients2D.push_back(coeff);
		}
}

void Gauss::Initialize(int ipc)
{
	m_IPC = ipc;
	calcCoefficients1D();
	calcCoefficients2D();
}
