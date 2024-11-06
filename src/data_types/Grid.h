#pragma once

#include <vector>

#include "Element.h"
#include "Node.h"

struct Grid
{
	std::vector<Element> Elements{};
	std::vector<Node> Nodes{};
};