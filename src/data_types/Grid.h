#pragma once

#include <vector>
#include <iostream>

#include "Element.h"
#include "Node.h"

struct Grid
{
	std::vector<Element> Elements{};
	std::vector<Node> Nodes{};

	std::vector<Node> getElementNodes(Element& element)
	{
		auto& nodeIDs = element.getNodesID();

		std::vector<Node> result(nodeIDs.size());
		for (size_t i = 0; i < nodeIDs.size(); ++i)
		{
			int nodeID = nodeIDs[i] - 1;
			result[i] = Nodes[nodeID];
		}
		return result;
	}
};