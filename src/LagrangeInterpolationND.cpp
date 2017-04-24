#include "../include/LagrangeInterpolationND.hpp"

LagrangeInterpolationND::~LagrangeInterpolationND()
{}

LagrangeInterpolationND::LagrangeInterpolationND(vector<int> sizes)
{
    m_d = sizes.size();
    m_points.resize(m_d);
    for (int i=0; i<m_d; i++)
        m_points[i].resize(sizes[i]);
}

void LagrangeInterpolationND::clear()
{
    m_path.clear();
    for (int i=0; i<m_d; i++)
        m_points[i].clear();
    m_points.clear();
    m_alphaTab.clear();
    m_curentNeighbours.clear();
}

/************************* Points **********************************************/
void LagrangeInterpolationND::showPoints()
{
    for (int i=0; i<m_d; i++)
    {
        cout << "Points d'interpolation suivant la direction " << i << ": ";
        for (int j=0; j<int(m_points[i].size()); j++)
            cout << m_points[i][j] << " ";
        cout << endl;
    }
}

/************************* Path ***********************************************/
void LagrangeInterpolationND::showPath()
{
    cout << "Chemin = ";
    set<IndiceND>::iterator it;
    for (it=m_path.begin(); it!=m_path.end(); it++)
    {
        cout << "(";
        for (int k=0; k<m_d; k++)
            cout << (*it)(k) << ",";
        cout << ")";
        if (it!=m_path.end())
            cout << " --> ";
    }
    cout << endl;
}
set<IndiceND>::iterator LagrangeInterpolationND::getLastIndiceInPath(vector<int> max)
{
    set<IndiceND>::iterator it = m_path.begin();
    bool found = false;
    while (!found)
    {
        found = true;
        for (int k=0; k<m_d; k++)
        {
            found = found && ((*it)(k)==max[k]);
            it++;
        }
    }
    return it;
}
bool LagrangeInterpolationND::indiceInPath(IndiceND& index)
{
    set<IndiceND>::iterator it = m_path.find(index);
    return (it!=m_path.end());
}
