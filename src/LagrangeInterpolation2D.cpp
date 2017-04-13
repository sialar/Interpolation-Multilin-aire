#include "../include/LagrangeInterpolation2D.hpp"

LagrangeInterpolation2D::~LagrangeInterpolation2D()
{}

LagrangeInterpolation2D::LagrangeInterpolation2D(int sizeX, int sizeY, int path)
{
    m_pointsX.resize(sizeX);
    m_pointsY.resize(sizeY);
    initPath(sizeX,sizeY,path);
    m_alphaTab.resize(sizeX);
    for (int i=0; i<sizeX; i++)
        m_alphaTab[i].resize(sizeY);
}

double LagrangeInterpolation2D::g(double x, double y)
{
    //return x*x*exp(y) + y*sin(2*x*x);
    return x*x + y*y*x + 2*y + 1;
}

void LagrangeInterpolation2D::initPath(int n, int m, int v /* 0 ou 1 ou 2*/)
{
    m_path.clear();
    indice2D index;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<((v%3)?m:i+1); j++)
        {
            index[0] = (v%3==1) ? i : j;
            index[1] = (v%3==0) ? i-j : ((v%3==1) ? j : i);
            m_path.push_back(index);
        }
    }
}

void LagrangeInterpolation2D::showPath()
{
    cout << "Chemin = ";
    for (int i=0; i<int(m_path.size()); i++)
    {
        cout << "(" << m_path[i][0] << "," << m_path[i][1] << ")";
        if (i<int(m_path.size())-1)
            cout << " --> ";
    }
    cout << endl;
}

void LagrangeInterpolation2D::showIndices()
{
    cout << "Indices = ";
    for (indice2D i : m_indices)
        cout << "(" << i[0] << "," << i[1] << ") ";
    cout << endl << endl;
}

void LagrangeInterpolation2D::updateIndices(int maxI, int maxJ)
{

    m_indices.clear();

    int _i = 0;
    indice2D i = m_path[_i];
    while (i[0]!=maxI || i[1]!=maxJ)
    {
        m_indices.push_back(i);
        _i++;
        i = m_path[_i];
    }
    /*
    cout << "Path = ";
    for (indice2D i : m_path)
        cout << "(" << i[0] << "," << i[1] << ") ";
    cout << endl << endl;

    cout << "Indices = ";
    for (indice2D i : m_indices)
        cout << "(" << i[0] << "," << i[1] << ") ";
    cout << endl << endl;
    */
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
    int _i = 0;
    while (m_path[_i][0]!=k1 || m_path[_i][1]!=k2)
    {
        lx = lagrangeBasisFunction_1D(m_path[_i][0],m_path[_i][0],x,0);
        ly = lagrangeBasisFunction_1D(m_path[_i][1],m_path[_i][1],y,1);
        sum += lx * ly * m_alphaTab[m_path[_i][0]][m_path[_i][1]];
        _i++;
    }
    return sum;
}

void LagrangeInterpolation2D::computeAllAlphaNu(int k1, int k2)
{
    m_alphaTab[0][0] = g(m_pointsX[0],m_pointsY[0]);
    int j1=0, j2=0;
    for (indice2D i : m_path)
    {
        // Ce qui suit, s'applique à chaque couple d'indices (l1,l2); l1 < k1, l2 < k2
        j1 = i[0]; j2 = i[1];
        updateIndices(j1,j2);
        m_alphaTab[j1][j2] = g(m_pointsX[j1],m_pointsY[j2]);
        for (indice2D l : m_indices)
        {
            m_alphaTab[j1][j2] -= m_alphaTab[l[0]][l[1]] * lagrangeBasisFunction_1D(l[0],l[0],m_pointsX[j1],0) *
                          lagrangeBasisFunction_1D(l[1],l[1],m_pointsY[j2],1);
        }
    }
}
