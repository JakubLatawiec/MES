#pragma once

#include <vector>
#include <unordered_map>
#include "../data_types/Node.h"

struct Coefficient1D { double X{}, W{}; };
struct Coefficient2D { Node Node{}; double SurfArea{}; };

class Gauss
{
private:
	static int m_IPC;
	static std::unordered_map<int, std::vector<Coefficient1D>> m_CoefficientsDatabase;

	static std::vector<Coefficient1D> m_Coefficients1D;
	static std::vector<Coefficient2D> m_Coefficients2D;

	static void calcCoefficients1D();
	static void calcCoefficients2D();

public:
	static void Initialize(int ipc);

	inline static const std::vector<Coefficient1D>& getIntegrationPoints1D() { return m_Coefficients1D; }
	inline static const std::vector<Coefficient2D>& getIntegrationPoints2D() { return m_Coefficients2D; }
};

