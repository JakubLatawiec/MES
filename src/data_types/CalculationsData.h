#pragma once

#include <stdexcept>
#include <vector>
#include "../utils/Gauss.h"

class CalculationsData
{
private:
	int m_IPC{};
	int m_IPC2D{};
	int m_SFC = 4;
	std::vector<Coefficient1D> m_IntegrationPoints{};
	std::vector<Coefficient2D> m_IntegrationPoints2D{};

	static CalculationsData m_Instance;
	bool m_IsInitialized = false;
	CalculationsData() = default;
	CalculationsData(const CalculationsData&) = delete;
	CalculationsData& operator=(const CalculationsData&) = delete;

public:
	static void Initialize(int ipc);
	static CalculationsData& getInstance();
	int getIPC() const;
	int getIPC2D() const;
	int getSFC() const;
	const std::vector<Coefficient1D>& getIntegrationPoints() const;
	const std::vector<Coefficient2D>& getIntegrationPoints2D() const;
};

