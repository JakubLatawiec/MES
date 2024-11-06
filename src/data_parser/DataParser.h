#pragma once

#include <vector>
#include <string>

#include "../data_types/GlobalData.h"
#include "../data_types/Grid.h"

class DataParser
{
private:
	static std::vector<std::string> split(const std::string& s, char delimeter);

public:
	static void ParseData(const std::string& filename, GlobalData& globalData, Grid& grid);
};

