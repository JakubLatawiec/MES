#include "Gauss.h"

#include <iostream>

std::unordered_map<int, std::vector<Coefficient1D>> Gauss::m_Coefficients =
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

std::vector<Coefficient1D> Gauss::GetIntegrationPoints1D(int ipc)
{
	return m_Coefficients[ipc];
}

std::vector<Coefficient2D> Gauss::GetIntegrationPoints2D(int ipc)
{
	std::vector<Coefficient2D> result{};
	const auto& points1D = m_Coefficients[ipc];

	for (const auto& csiPoint : points1D)
		for (const auto& etaPoint : points1D)
		{
			Coefficient2D coeff;
			coeff.Node.csi = csiPoint.X;
			coeff.Node.eta = etaPoint.X;
			coeff.SurfArea = csiPoint.W * etaPoint.W;
			result.push_back(coeff);
		}

	return result;
}
