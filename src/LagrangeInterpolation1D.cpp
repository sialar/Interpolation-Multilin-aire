#include "../include/LagrangeInterpolation1D.hpp"

LagrangeInterpolation1D::~LagrangeInterpolation1D()
{}

LagrangeInterpolation1D::LagrangeInterpolation1D(int size)
{
    m_points.resize(size);
    allLiMinus1Fi.resize(size+1);
}

double LagrangeInterpolation1D::g(double y)
{
    return y*y*y + y*y + y + 1;
}

double LagrangeInterpolation1D::lagrangeBasisFunction_1D(int j, int k, double y)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        if (i!=j)
            prod *= ( (y-m_points[i]) / (m_points[j]-m_points[i]) );
    return prod;
}

// 1 ère methode
double LagrangeInterpolation1D::lagrangeInterpolation_1D_simple(double y)
{
    double sum = 0;
    for (int i=0; i<int(m_points.size()); ++i)
        sum += g(m_points[i]) * lagrangeBasisFunction_1D(i,m_points.size(),y);
    return sum;
}

// 2 ère methode (version récurssive ~ trés couteuse : calcul des même quantités plusieurs fois)
double LagrangeInterpolation1D::lagrangeInterpolation_1D_recursive(double y, int k)
{
    if (k==-1) return 0;
    double sum = 0;
    for (int i=0; i<=k; ++i)
        sum += lagrangeBasisFunction_1D(i,i,y) * (g(m_points[i]) - lagrangeInterpolation_1D_recursive(m_points[i],i-1));
    return sum;
}

// 2 ère methode (version itérative couteuse en terme de mèmoire)
double LagrangeInterpolation1D::lagrangeInterpolation_1D_iterative(double y,int k)
{
    double sum = 0;
    for (int i=0; i<k; i++)
        sum += lagrangeBasisFunction_1D(i,i,y) * (g(m_points[i]) - allLiMinus1Fi[i]);
    return sum;
}

void LagrangeInterpolation1D::computeLiMinus1Fi(int k)
{
    // Calcul de tous les I_{l} f(t_{l+1}) et stockage dans le tableau allLiMinus1Fi
    allLiMinus1Fi[0] = 0; //I_{-1} f(t_0)
    for (int l=1; l<=k; ++l)
    {
        allLiMinus1Fi[l] = 0;
        for (int i=0; i<l; ++i)   // I_{l-1} f(t_{l}) = sum_{i=0}^{l-1} ( (g(t_i)-I{i-1} f(t_i)) * h_i(t_{l})
            allLiMinus1Fi[l] += (g(m_points[i]) - allLiMinus1Fi[i]) *  lagrangeBasisFunction_1D(i,i,m_points[l]);
    }
}
