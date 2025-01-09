#include "CalculationsData.h"

CalculationsData CalculationsData::m_Instance;

void CalculationsData::Initialize(int ipc)
{
    if (m_Instance.m_IsInitialized)
        throw std::runtime_error("CalculationsData is already initialized!");

    if (ipc <= 0)
        throw std::invalid_argument("IPC must be greater than 0!");

    m_Instance.m_IPC = ipc;
    m_Instance.m_IPC2D = static_cast<int>(pow(ipc, 2));
    m_Instance.m_IntegrationPoints = Gauss::GetIntegrationPoints1D(ipc);
    m_Instance.m_IntegrationPoints2D = Gauss::GetIntegrationPoints2D(ipc);
    m_Instance.m_IsInitialized = true;
}

CalculationsData& CalculationsData::getInstance()
{
    if (!m_Instance.m_IsInitialized)
        throw std::runtime_error("CalculationsData has not been initialized!");

    return m_Instance;
}

int CalculationsData::getIPC() const
{
    return m_IPC;
}

int CalculationsData::getIPC2D() const
{
    return m_IPC2D;
}

int CalculationsData::getSFC() const
{
    return m_SFC;
}

const std::vector<Coefficient1D>& CalculationsData::getIntegrationPoints() const
{
    return m_IntegrationPoints;
}

const std::vector<Coefficient2D>& CalculationsData::getIntegrationPoints2D() const
{
    return m_IntegrationPoints2D;
}
