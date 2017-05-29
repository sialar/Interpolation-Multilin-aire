#include "../include/Interpolation.hpp"

Interpolation::Interpolation(int d, int nIter) : m_d(d), m_maxIteration(nIter)
{
    m_interpolationPoints.resize(m_d);
}

void Interpolation::displayInterpolationPointsInEachDirection()
{
    vector<double>::iterator it;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_interpolationPoints[i].size() << " points in direction " << i << " : { ";
        for (it=m_interpolationPoints[i].begin(); it!=m_interpolationPoints[i].end(); it++)
            cout << *it << " ";
        cout << "}" << endl;
    }
}

void Interpolation::displayInterpolationMultiVariatePoints()
{
    cout << " - Interpolation nodes: { ";
    for (MultiVariatePoint<double> x : m_interpolationNodes)
        cout << x << " ";
    cout << "}" << endl;
}
