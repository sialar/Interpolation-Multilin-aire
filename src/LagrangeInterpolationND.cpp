#include "../include/LagrangeInterpolationND.hpp"

LagrangeInterpolationND::~LagrangeInterpolationND()
{}

LagrangeInterpolationND::LagrangeInterpolationND(int size)
{
    m_points.resize(size);
}
