#include "../include/LagrangeInterpolation2D.hpp"

LagrangeInterpolation2D::~LagrangeInterpolation2D()
{}

LagrangeInterpolation2D::LagrangeInterpolation2D(int sizeX, int sizeY, int path)
{
    m_pointsX.resize(sizeX);
    m_pointsY.resize(sizeY);
    initPath(sizeX,sizeY,path);
    showPath();
    values.resize(sizeX);
    for (int i=0; i<sizeX; i++)
        values[i].resize(sizeY);
}

double LagrangeInterpolation2D::g(double x, double y)
{
    return x*x*exp(y) + y*sin(2*x*x);
}

void LagrangeInterpolation2D::initPath(int n, int m, int v /* 0 ou 1 */)
{

    // path version 0
    indice2D index;
    int epsilon = 0;
    int i = 0;
    for (int j=0; j<((v%2)?m:n); j++)
    {
        if (epsilon%2==0)
        {
            i = 0;
            while (i<((v%2)?n:m))
            {
                index[0] = ((v%2)?i:j);
                index[1] = ((v%2)?j:i);
                m_path.push_back(index);
                i++;
            }
            epsilon++;
        }
        else
        {
            i =  ((v%2)?n:m)-1;
            while (i>=0)
            {
                index[0] = ((v%2)?i:j);
                index[1] = ((v%2)?j:i);
                m_path.push_back(index);
                i--;
            }
            epsilon++;
        }
    }
    //m_path.erase(m_path.begin());
    if (v==0)
    {
        // A verifier (il existe des sauts)
        m_path.clear();
        indice2D index;
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<m; j++)
            {
                index[0]=i; index[1]=j;
                m_path.push_back(index);
            }
        }
    }
}

void LagrangeInterpolation2D::showPath()
{
    cout << "Path = ";
    for (indice2D i : m_path)
        cout << "(" << i[0] << "," << i[1] << ") ";
    cout << endl << endl;
}

void LagrangeInterpolation2D::showIndices()
{
    cout << "Indices = ";
    for (indice2D i : m_indices)
        cout << "(" << i[0] << "," << i[1] << ") ";
    cout << endl << endl;
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
        sum += lx * ly * (g(m_pointsX[m_path[_i][0]],m_pointsY[m_path[_i][1]]) - values[m_path[_i][0]][m_path[_i][1]]);
        _i++;
    }
    return sum;
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

void LagrangeInterpolation2D::computeValues(int k1, int k2)
{
    values[0][0] = 0;
    int j1=0, j2=0;
    for (indice2D i : m_path)
    {
        // Ce qui suit, s'applique à chaque couple d'indices (l1,l2); l1 < k1, l2 < k2
        j1 = i[0]; j2 = i[1];
        updateIndices(j1,j2);
        //showIndices();
        for (indice2D l : m_indices)
        {
            values[j1][j2] += (g(m_pointsX[l[0]],m_pointsY[l[1]]) - values[l[0]][l[1]]) *
                lagrangeBasisFunction_1D(l[0],l[0],m_pointsX[j1],0) * lagrangeBasisFunction_1D(l[1],l[1],m_pointsY[j2],1);
        }
    }
}
