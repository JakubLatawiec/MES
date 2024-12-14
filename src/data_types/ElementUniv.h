#pragma once

#include "../utils/Matrix.h"
#include "../utils/Gauss.h"
#include <array>
#include <iostream>


struct ElementUniv
{
	Matrix CsiDerivative;
	Matrix EtaDerivative;
	Matrix PcN{};

	void calcPcN(int ipc);

	//TODO:
	//SEPARATE ELEMENT UNIV CALC FUNCTIONS
};