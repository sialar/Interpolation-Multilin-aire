#include "../include/LagrangeInterpolation2D.hpp"

LagrangeInterpolation2D::~LagrangeInterpolation2D()
{}

LagrangeInterpolation2D::LagrangeInterpolation2D(int sizeX, int sizeY, int path)
{
    m_pointsX.resize(sizeX);
    m_pointsY.resize(sizeY);
    choosePath(sizeX,sizeY,path);
    m_alphaTab.resize(sizeX);
    for (int i=0; i<sizeX; i++)
        m_alphaTab[i].resize(sizeY);
}

void LagrangeInterpolation2D::clear()
{
    m_pointsX.clear();
    m_pointsY.clear();
    m_path.clear();
    for (int i=0; i<int(m_pointsX.size()); i++)
        m_alphaTab[i].clear();
    m_alphaTab.clear();
    m_curSetInAIAlgo.clear();
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

void LagrangeInterpolation2D::showAlphaTab()
{
    cout << "AlphaTab = " << endl;
    for (int i=0; i<int(m_alphaTab.size()); i++)
    {
        for (int j=0; j<int(m_alphaTab[0].size()); j++)
            cout << m_alphaTab[i][j] << " ";
        cout << endl;
    }
}

double LagrangeInterpolation2D::g(double x, double y)
{
    return sin(sqrt(x*x + y*y));
    //return x*x*exp(y) + pow(y*sin(2*x*x),4);
    //return x*x + y*y*x + 2*y + 1;
}

void LagrangeInterpolation2D::choosePath(int n, int m, int v /* 0 ou 1 ou 2*/)
{
    m_path.clear();
    if (v >= 0)
    {
        indice2D index;
        for (int i=0; i<n; i++)
            for (int j=0; j<((v%3)?m:i+1); j++)
            {
                index =  { (v%3==1) ? i : j , (v%3==0) ? i-j : ((v%3==1) ? j : i)};
                m_path.push_back(index);
            }
    }
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
        for (int j=0; j<int(m_pointsY.size()); j++)
        {
            lx = lagrangeBasisFunction_1D(i,m_pointsX.size(),x,0);
            ly = lagrangeBasisFunction_1D(j,m_pointsY.size(),y,1);
            z += g(x,y)*lx*ly;
        }
    return z;
}

double LagrangeInterpolation2D::lagrangeInterpolation_2D_iterative(double x, double y)
{
    double sum = 0, lx = 0, ly = 0;
    for (int i=0; i<int(m_path.size()); i++)
    {
        lx = lagrangeBasisFunction_1D(m_path[i][0],m_path[i][0],x,0);
        ly = lagrangeBasisFunction_1D(m_path[i][1],m_path[i][1],y,1);
        sum += lx * ly * m_alphaTab[m_path[i][0]][m_path[i][1]];
    }
    return sum;
}

int LagrangeInterpolation2D::getIndiceInPath(int maxI, int maxJ)
{
    int i = 0;
    while (m_path[i][0]!=maxI || m_path[i][1]!=maxJ)
        i++;
    return i;
}

void LagrangeInterpolation2D::computeAllAlphaNu()
{
    for (indice2D i : m_path)
        computeOneAlphaNu(i);
}

double LagrangeInterpolation2D::computeOneAlphaNu(indice2D nu)
{
    int j1 = nu[0], j2 = nu[1];
    int i = getIndiceInPath(j1,j2);
    m_alphaTab[j1][j2] = g(m_pointsX[j1],m_pointsY[j2]);
    for (int k=0; k<i; k++)
    {
        indice2D l = m_path[k];
        m_alphaTab[j1][j2] -= m_alphaTab[l[0]][l[1]] * lagrangeBasisFunction_1D(l[0],l[0],m_pointsX[j1],0) *
                      lagrangeBasisFunction_1D(l[1],l[1],m_pointsY[j2],1);
    }
    return m_alphaTab[j1][j2];
}

double LagrangeInterpolation2D::computeLastAlphaNu(indice2D nu)
{
    int j1 = nu[0], j2 = nu[1];
    m_alphaTab[j1][j2] = g(m_pointsX[j1],m_pointsY[j2]);
    for (int k=0; k<int(m_path.size()); k++)
    {
        indice2D l = m_path[k];
        m_alphaTab[j1][j2] -= m_alphaTab[l[0]][l[1]] * lagrangeBasisFunction_1D(l[0],l[0],m_pointsX[j1],0) *
                      lagrangeBasisFunction_1D(l[1],l[1],m_pointsY[j2],1);
    }
    return m_alphaTab[j1][j2];
}

vector<indice2D> LagrangeInterpolation2D::getCurentNeighbours()
{
    vector<indice2D> neighbours;
    indice2D index;
    for (indice2D nu : m_curSetInAIAlgo)
        if (!indexInPath(nu))
            for (int i=0; i<2; i++)
                if (nu[i]>0)
                {
                    index[i] = nu[i]-1;
                    index[(i+1)%2] = nu[(i+1)%2];
                    if (indexInPath(index))
                        neighbours.push_back(nu);
                }
    return neighbours;
}

/*
Ici on met à jour la liste d'indices courante en fonction de ce qu'on a dans le vecteur path.
Avec l'algo AI, on veut en meme temps construire le chemin et calculer les Alpha
On ne va plus mettre a jour la liste d'indices a partir du chemin puisqu'il n'est pas totalement construit.
On va donc recupere l'information sur les indices qui interviennent dans le calcul du alpha courant grace à l'algo AI
En effet le nouveau point ajouté ne sera pas le courant dans le chemin comme dans la methode precedante mais il sera le point renvoyé par l'algo AI.
*/
void LagrangeInterpolation2D::buildPathWithAIAlgo(int n, int m, bool isLejaSeq)
{
    m_curSetInAIAlgo.clear();
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            m_curSetInAIAlgo.push_back({i,j});
    m_curSetInAIAlgo.erase (m_curSetInAIAlgo.begin());

    indice2D index = {0,0};
    m_path.clear();
    m_path.push_back(index);
    vector<indice2D> neighbours;

    double val, max;
    indice2D argmax;

    while (int(m_path.size()) < n*m)
    {
        //cout << "size of path = " << m_path.size() << "\n";
        max = -numeric_limits<double>::max();
        neighbours = getCurentNeighbours();

        // Test the funcion getCurentNeighbours() --> OK !!
        //cout << "Neighbours of {";
        //displayVector(m_path);
        //cout << "} are {";
        //displayVector(neighbours);
        //cout << "}" << endl;
        //////////////////////////////////////////////////

        for (indice2D nu : neighbours)
        {
            if (isLejaSeq)
            {
                // Note cette fonction permet de calculer alphaI sans utiliser le path
                val = computeLastAlphaNu(nu);
                //cout << val << endl;
            }
            else
            {
                //val = ..
            }
            if (val >= max)
            {
                max = val;
                argmax = nu;
            }
        }
        //cout << "new indice chosen = (" << argmax[0] << "," << argmax[1] << ")" << m_path.size() << "\n";
        m_path.push_back(argmax);

    }
}

bool LagrangeInterpolation2D::indexInPath(indice2D index)
{
    for (indice2D _i : m_path)
        if (_i[0]==index[0] && _i[1]==index[1])
            return true;
    return false;
}
