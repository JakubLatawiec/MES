#pragma once

struct Node 
{
    union
    {
        struct { double x, y; };
        struct { double csi, eta; };
    };

    bool isBorderCondition = false;

    Node() : x(0.0), y(0.0) {}
    Node(double x, double y) : x(x), y(y) {}
};