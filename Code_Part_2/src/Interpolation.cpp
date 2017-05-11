#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{}

Interpolation::Interpolation(int size)
{
    m_alphaTab.resize(size+1);

    m_tree = new BinaryTree();
    m_tree->initTree(int(log(size)/log(2))-1);
    //m_tree->displayBinaryTree();

    m_points.push_back(-1);
    m_points.push_back(1);
    m_points.push_back(0);
    int e = 0, i = 3;
    double denom;
    while (i<size)
    {
        e++;
        denom = pow(2,e);
        for (int j=denom-1; j>=0; j--)
            if (j%2)
                m_points.push_back(-j/denom);
        for (int j=0; j<denom; j++)
            if (j%2)
                m_points.push_back(j/denom);
        i += denom;
    }
    for (double d : m_points)
        cout << d << " ";
    cout << endl;
    /*
    // Test creation des points
    separateur();
    double val;
    int p=0;
    while (p<size)
    {
        val = m_points[p];
        cout << "Indice found: " << getIndice(val) << endl;
        cout << "Checking ...\t" <<  m_points[getIndice(val)] << " should be " << val << endl << endl;
        p++;
    }
    */
}

double Interpolation::piecewiseFunction_1D(int k, double t)
{
    double sup, inf, tk = m_points[k];
    m_tree->searchNode(tk,&sup,&inf);
    if (t <= tk && t >= inf) return ((t-inf)/(tk-inf));
    else if (t >= tk && t <= sup) return ((t-sup)/(tk-sup));
    else return 0;
}

bool Interpolation::indiceInPath(int index)
{
    for (int i : m_path)
        if (index==i)
            return true;
    return false;
}

void Interpolation::updateNextPoints(int i)
{
    m_nextPoints.remove(i);
    if (i == 0 && !indiceInPath(1)) m_nextPoints.push_back(1);
    else if (i == 1 && !indiceInPath(2)) m_nextPoints.push_back(2);
    else
    {
        double inf, sup;
        Node* node = m_tree->searchNode(m_points[i],&inf,&sup);
        int inf_index, sup_index;
        if (node->left())
        {
            inf = node->left()->key();
            inf_index = getIndice(inf);
            if (!indiceInPath(inf_index))
                m_nextPoints.push_back(getIndice(inf));
        }
        if (node->right())
        {
            sup = node->right()->key();
            sup_index = getIndice(sup);
            if (!indiceInPath(sup_index))
                m_nextPoints.push_back(getIndice(sup));
        }
    }
}

double Interpolation::interpolation_iterative(double y, int k, bool debug)
{
    double sum = 0;
    m_path.push_back(0);
    m_alphaTab[0] = g(m_path[0]);
    double val, max;
    int argmax = 0, iteration = 0;
    while (iteration < k)
    {
        max = -numeric_limits<double>::max();
        updateNextPoints(argmax);
        if (debug)
        {
          displayPath();
          displayNextPoints();
        }
        for (int p : m_nextPoints)
        {
            val =  computeLastAlphaI(p);
            if (debug) cout << "(" << p << "," << m_points[p] << ") : " << val << " | ";
            if (val >= max)
            {
                max = val;
                argmax = p;
            }
        }
        if (debug) cout << endl << "(" << argmax << "," << m_points[argmax] << ")" << endl << endl;
        m_path.push_back(argmax);
        iteration++;
    }
    for (int i : m_path)
        sum += piecewiseFunction_1D(i,y) * m_alphaTab[i];
    return sum;
}

int Interpolation::getIndice(double l)
{
    if (l == -1) return 0;
    if (l == 1)  return 1;
    if (l == 0)  return 2;
    double temp = l;
    int n = 0;
    while (floor(temp) != temp)
    {
        temp *= 2;
        n++;
    }
    int start = 2;
    for (int k=1; k<n; k++)
        start += pow(2,n-k);
    double a = 0.5 * (temp + pow(2,n)) + 1;
    return start + int(a);
}

double Interpolation::computeLastAlphaI(int i)
{
    double res = g(m_points[i]);
    for (double l : m_path)
        res -= m_alphaTab[l] * piecewiseFunction_1D(l, m_points[i]);
    m_alphaTab[i] = res;
    return res;
}

void Interpolation::displayPath()
{
  cout << "Chemin = ";
  for (int i=0; i<int(m_path.size()); i++)
  {
        cout << "(" << m_path[i] << "," << m_points[m_path[i]] << ")";
    if (i != int(m_path.size()-1))
        cout << " --> ";
  }
  cout << endl;
}

void Interpolation::displayAlphaTab()
{
    cout << "Values of alpha (" << m_alphaTab.size() << ") :" << endl;
    for (int i : m_path)
        cout << i << ": " << m_alphaTab[i] << endl;
}

void Interpolation::displayNextPoints()
{
    cout << "Next points (" << m_nextPoints.size() << ") = ";
    for (int i: m_nextPoints)
    cout << "(" << i << "," << m_points[i] << ") ";
    cout << endl;
}
