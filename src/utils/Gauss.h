#pragma once

#include <vector>

#include "../data_types/Node.h"

class Gauss
{
private:
	struct Coefficient1D { double X, W; };
	struct Coefficient2D { Node Node; double SurfArea; };

public:
	static std::vector<Coefficient1D> GetIntegrationPoints1D(int npc);
	static std::vector<Coefficient2D> GetIntegrationPoints2D(int npc);
};

