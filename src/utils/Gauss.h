#pragma once

#include <vector>
#include <unordered_map>

#include "../data_types/Node.h"

struct Coefficient1D 
{
	double X{}, W{};
	Coefficient1D() = default;
};
struct Coefficient2D { Node Node{}; double SurfArea{}; };

class Gauss
{
private:
	static std::unordered_map<int, std::vector<Coefficient1D>> m_Coefficients;
public:
	static std::vector<Coefficient1D> GetIntegrationPoints1D(int npc);
	static std::vector<Coefficient2D> GetIntegrationPoints2D(int npc);
};

