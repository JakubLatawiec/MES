#include "Gauss.h"

//TODO: Wyliczanie analityczne
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
			{-0.861136311594053, 0.347854845137454},
			{-0.339981043584856, 0.652145154862546},
			{ 0.339981043584856, 0.652145154862546},
			{ 0.861136311594053, 0.347854845137454}
		}
	},
	{5, {
			{-0.906179845938664, 0.236926885056189},
			{-0.538469310105683, 0.478628670499366},
			{0.0,                0.568888888888889},
			{ 0.538469310105683, 0.478628670499366},
			{ 0.906179845938664, 0.236926885056189}
		}
	}
};

std::vector<Coefficient1D> Gauss::GetIntegrationPoints1D(int npc)
{
	std::vector<Coefficient1D> res{};

	switch (npc)
	{
	case 1:
		res =
		{
			{0.0, 2.0}
		};
		break;
	case 2:
		res =
		{
			{-1.0 / sqrt(3.0), 1.0},
			{1.0 / sqrt(3.0), 1.0}
		};
		break;
	case 3:
		res =
		{
			{-sqrt(3.0 / 5.0), 5.0 / 9.0},
			{0.0, 8.0 / 9.0},
			{sqrt(3.0 / 5.0), 5.0 / 9.0}
		};
		break;
	case 4:
		res =
		{
			{-0.861136, 0.347855},
			{-0.339981, 0.652145},
			{0.339981, 0.652145},
			{0.861136, 0.347855}
		};
		break;
	default:
		break;
	}

	return res;
}

std::vector<Coefficient2D> Gauss::GetIntegrationPoints2D(int npc)
{
	std::vector<Coefficient2D> res{};

	if (npc == 4)
	{
		double w = 1.0;
		double a = 1.0 / sqrt(3.0);
		res =
		{
			{Node(-a, -a), w},
			{Node(a, -a), w},
			{Node(-a, a), w},
			{Node(a, a), w}
		};
	}
	else if (npc == 9)
	{
		double w1 = 5.0 / 9.0;
		double w2 = 8.0 / 9.0;
		double a = sqrt(3.0 / 5.0);

		res =
		{
			{Node(-a, -a), w1 * w1},
			{Node(a, -a), w1 * w1},
			{Node(-a, a), w1 * w1},
			{Node(a, a), w1 * w1},
			{Node(0, -a), w2 * w1},
			{Node(0, a), w2 * w1},
			{Node(-a, 0), w1 * w2},
			{Node(a, 0), w1 * w2},
			{Node(0, 0), w2 * w2},
		};
	}

	return res;
}
