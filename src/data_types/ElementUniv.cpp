#include "ElementUniv.h"
#include "../utils/Gauss.h"

void ElementUniv::calcPcN(int ipc)
{
	PcN = Matrix(ipc, 4);
	auto integrationPoints = Gauss::GetIntegrationPoints2D(ipc);

	for (size_t i = 0; i < ipc; ++i)
	{
		double csi = integrationPoints[i].Node.csi;
		double eta = integrationPoints[i].Node.eta;

		PcN(i, 0) = 0.25 * (1 - csi) * (1 - eta);
		PcN(i, 1) = 0.25 * (1 + csi) * (1 - eta);
		PcN(i, 2) = 0.25 * (1 + csi) * (1 + eta);
		PcN(i, 3) = 0.25 * (1 - csi) * (1 + eta);
	}
}
