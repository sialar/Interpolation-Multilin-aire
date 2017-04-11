#include "../include/LagrangeInterpolation2D.hpp"

LagrangeInterpolation2D::~LagrangeInterpolation2D()
{}

LagrangeInterpolation2D::LagrangeInterpolation2D(int sizeX, int sizeY)
{
    m_pointsX.resize(sizeX);
    m_pointsY.resize(sizeY);
    allLiMoins1Fi.resize(sizeX);
    for (int i=0; i<sizeX; i++)
        allLiMoins1Fi[i].resize(sizeY);
}

double LagrangeInterpolation2D::g(double x, double y)
{
    return x*x + y*y;
}

double LagrangeInterpolation2D::lagrangeBasisFunction_1D(int j, int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        if (i!=j)
            prod *= (axis==0) ? ( (t-m_pointsX[i]) / (m_pointsX[j]-m_pointsX[i]) )
                : ( (t-m_pointsY[i]) / (m_pointsY[j]-m_pointsY[i]) );
    return prod;
}

double LagrangeInterpolation2D::lagrangeInterpolation_2D_simple(double x, double y)
{
    double z = 0, lx = 0, ly = 0;
    for (int i=0; i<int(m_pointsX.size()); i++)
    {
        for (int j=0; j<int(m_pointsY.size()); j++)
        {
            lx = lagrangeBasisFunction_1D(i,m_pointsX.size(),x,0);
            ly = lagrangeBasisFunction_1D(j,m_pointsY.size(),y,1);
            z += g(x,y)*lx*ly;
        }
    }
    return z;
}

double LagrangeInterpolation2D::lagrangeInterpolation_2D_iterative(double x, double y, int k1, int k2)
{
    double sum = 0, lx = 0, ly = 0;
    for (int i=0; i<k1; i++)
    {
        for (int j=0; j<k2; j++)
        {
            lx = lagrangeBasisFunction_1D(i,i,x,0);
            ly = lagrangeBasisFunction_1D(j,j,y,1);
            sum += lx * ly * (g(m_pointsX[i],m_pointsY[j]) - allLiMoins1Fi[i][j]);
        }
    }
    return sum;
}

void LagrangeInterpolation2D::fillIndicesWithoutMax(int maxI, int maxJ)
{
    m_indices.clear();
    indice2D index;
    for (int i=0; i<maxI; i++)
    {
        for (int j=0; j<=maxJ; j++)
        {
            index[0]=i; index[1]=j;
            m_indices.push_back(index);
        }
    }
    for (int j=0; j<maxJ; j++)
    {
        index[0]=maxI; index[1]=j;
        m_indices.push_back(index);
    }
}

void LagrangeInterpolation2D::computeLiMinus1Fi(int k1, int k2)
{
    allLiMoins1Fi[0][0] = 0;
    for (int j1=0; j1<k1; ++j1)
    {
        for (int j2=0; j2<k2; ++j2)
        {
            // Ce qui suit, s'applique à chaque couple d'indices (l1,l2); l1 < k1, l2 < k2
            fillIndicesWithoutMax(j1,j2);
            for (indice2D l : m_indices)
            {
                allLiMoins1Fi[j1][j2] += (g(m_pointsX[l[0]],m_pointsY[l[1]]) - allLiMoins1Fi[l[0]][l[1]]) *
                    lagrangeBasisFunction_1D(l[0],l[0],m_pointsX[j1],0) * lagrangeBasisFunction_1D(l[1],l[1],m_pointsY[j2],1);
            }
        }
    }
}
