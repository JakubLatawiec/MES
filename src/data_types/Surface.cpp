#include "Surface.h"

#include "../utils/Gauss.h"
#include <iostream>

void Surface::calcSurface(int npc)
{
	auto integrationPoints = Gauss::getIntegrationPoints1D();

	for (int side = 0; side < 4; ++side)
	{
		double csi = 0;
		double eta = 0;

		this->HbcMatrixes[side] = Matrix(4, 4);
		this->PVectors[side] = Matrix(4, 1);

		Matrix pcN = Matrix(npc, 4);

		for (int i = 0; i < npc; ++i)
		{
			switch (side)
			{
			case 0:
				csi = integrationPoints[i].X;
				eta = -1;
				break;
			case 1:
				csi = 1;
				eta = integrationPoints[i].X;
				break;
			case 2:
				csi = integrationPoints[i].X;
				eta = 1;
				break;
			case 3:
				csi = -1;
				eta = integrationPoints[i].X;
				break;
			}

			pcN(i, 0) = 0.25 * (1 - csi) * (1 - eta);
			pcN(i, 1) = 0.25 * (1 + csi) * (1 - eta);
			pcN(i, 2) = 0.25 * (1 + csi) * (1 + eta);
			pcN(i, 3) = 0.25 * (1 - csi) * (1 + eta);
		}

		for (int i = 0; i < npc; ++i)
		{
			this->HbcMatrixes[side] = this->HbcMatrixes[side] + (pcN.getRow(i).Transpose() * pcN.getRow(i) * integrationPoints[i].W);
			this->PVectors[side] = this->PVectors[side] + (pcN.getRow(i).Transpose() * integrationPoints[i].W);
		}
			
			
	}
}
