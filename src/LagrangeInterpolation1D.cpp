#include "../include/LagrangeInterpolation1D.hpp"

LagrangeInterpolation1D::~LagrangeInterpolation1D()
{}

LagrangeInterpolation1D::LagrangeInterpolation1D(int size)
{
    m_points.resize(size);
    m_alphaTab.resize(size+1);
}

double LagrangeInterpolation1D::g(double y)
{
    return y*y;
    //return y*y*y + pow((y+1)*sin(y),3) + y*exp(2*y) + 1;
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

// 1 ère methode : méthode simple
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
        sum += lagrangeBasisFunction_1D(i,i,y) * m_alphaTab[i];
    return sum;
}

void LagrangeInterpolation1D::computeAllAlphaI(int k)
{
    // Calcul de tous les alpha_{j} et stockage dans le tableau m_alphaTab
    m_alphaTab[0] = g(m_points[0]); //alpha_0
    for (int j=1; j<=k; ++j)
    {
        m_alphaTab[j] = g(m_points[j]);
        for (int l=0; l<j; ++l)   // alpha_j = g(t_j) - sum_{l=0}^{j-1} ( alpha_l) * h_l(t_j)
            m_alphaTab[j] -= m_alphaTab[l] *  lagrangeBasisFunction_1D(l,l,m_points[j]);
    }
}
