#include "../include/PiecewiseInterpolation.hpp"

PiecewiseInterpolation::PiecewiseInterpolation(int d, int nIter, int method) :
    Interpolation(d, nIter), m_method(method)
{
    m_trees.resize(m_d);
    for (int i=0; i<m_d; i++)
      m_trees[i] = make_shared<BinaryTree>();
  }

/************************* Data Points ****************************************/
MultiVariatePoint<double> PiecewiseInterpolation::getPoint(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<double> point(m_d,0.0);
    for (int i=0; i<m_d; i++)
        point(i) = BinaryTree::getValueFromCode((*nu)(i));
    return point;
}
void PiecewiseInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
{
    m_interpolationNodes.push_back(p);

    bool found;
    for (int i=0; i<m_d; i++)
    {
        found = false;
        for (double x : m_interpolationPoints[i])
            if (x == p(i)) found = true;
        if (!found)
        {
            m_interpolationPoints[i].push_back(p(i));
            m_trees[i]->addNode(p(i));
        }
    }
}
void PiecewiseInterpolation::computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis)
{
    vector<double> higherPoints, lowerPoints;
    for (double x : m_interpolationPoints[axis])
    {
        if (x>t) higherPoints.push_back(x);
        else if (x<t) lowerPoints.push_back(x);
        else break;
        // break when t is found, do not look at children when looking for boundaries
        // because children were constructed after t
    }
    if (higherPoints.size())
        *sup = *min_element(higherPoints.begin(),higherPoints.end());
    else *sup = t;
    if (lowerPoints.size())
        *inf = *max_element(lowerPoints.begin(),lowerPoints.end());
    else *inf = t;

    // Or using trees
    m_trees[axis]->searchNode(t,inf,sup,false);
}

/******************************************************************************/


/************************* AI algo ********************************************/
bool alphaLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return abs(nu->getAlpha()) < abs(mu->getAlpha());
}
bool ageLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return nu->getWaitingTime() < mu->getWaitingTime();
}
MultiVariatePointPtr<string> PiecewiseInterpolation::maxElement(int iteration)
{
    return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),(iteration%4) ? alphaLess : ageLess);
}
MultiVariatePointPtr<string> PiecewiseInterpolation::getFirstMultivariatePoint()
{
    MultiVariatePointPtr<string> nu = make_shared<MultiVariatePoint<string>>(m_d,"");
    return nu;
}
void PiecewiseInterpolation::updateCurentNeighbours(MultiVariatePointPtr<string> nu)
{
  m_curentNeighbours.remove(nu);
  for (MultiVariatePointPtr<string> mu : m_curentNeighbours)
      mu->incrWaitingTime();

  vector<string> childrenCodes;
  for (int i=0; i<m_d; i++)
  {
      childrenCodes = BinaryTree::computeChildrenCodes((*nu)(i));
      for (int k=0; k<int(childrenCodes.size()); k++)
      {
          MultiVariatePointPtr<string> mu = make_shared<MultiVariatePoint<string>>(*nu);
          (*mu)(i) = childrenCodes[k];
          if (isCorrectNeighbourToCurentPath(mu))
              m_curentNeighbours.push_back(mu);
      }
  }
}
bool PiecewiseInterpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu)
{
    bool isNeighbour[m_d];
    for (int i=0; i<m_d; i++)
        isNeighbour[i] = true;

    MultiVariatePoint<string> mu(*nu);
    string temp;
    for (int i=0; i<m_d; i++)
        if ((*nu)(i).compare("")!=0)
        {
            temp = mu(i);
            mu(i) = BinaryTree::getParentCode(temp);
            isNeighbour[i] = indiceInPath(mu) && !indiceInNeighborhood(mu);
            mu(i) = temp;
        }
    bool res = true;
    for (int i=0; i<m_d; i++)
        res = res && isNeighbour[i];
    return res;
}
bool PiecewiseInterpolation::indiceInNeighborhood(MultiVariatePoint<string> index)
{
    for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
        if (Utils::equals(index,*nu))
            return true;
    return false;
}
bool PiecewiseInterpolation::indiceInPath(MultiVariatePoint<string> index)
{
    for (MultiVariatePointPtr<string> nu : m_path)
    if (Utils::equals(index,*nu))
        return true;
    return false;
}
/******************************************************************************/

/************************* Display functions **********************************/
void PiecewiseInterpolation::displayTrees()
{
    for (int i=0; i<m_d; i++)
    {
        cout << " - Binary Tree in direction " << i << " :" << endl;
        m_trees[i]->displayBinaryTree();
        cout << endl;
    }
}
void PiecewiseInterpolation::storeInterpolationBasisFunctions()
{
  ofstream file("data/basis_functions.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==1)
    {
      vector<double> x;
      int nbPoints = int(m_testPoints.size());
      for (int i=0; i<nbPoints; i++)
          x.push_back(m_testPoints[i](0));
      sort(x.begin(), x.end());
      file << m_path.size() << " " << nbPoints << " " << m_method << endl;
      MultiVariatePoint<double> p;
      for (int j=0; j<nbPoints; j++)
      {
        file << x[j];
        for (MultiVariatePointPtr<string> nu : m_path)
          file << " " <<  basisFunction_1D((*nu)(0),x[j],0);
        p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
        file << " " <<  Utils::gNd(p);
        file << endl;
      }
      for (MultiVariatePointPtr<string> nu : m_path)
      file << nu->getAlpha() << " ";
      file << endl;
      for (MultiVariatePoint<double> nu : m_interpolationNodes)
      file << nu(0) << " ";
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
/******************************************************************************/

/*********************** Interpolation ****************************************/
double PiecewiseInterpolation::basisFunction_1D(string code, double t, int axis)
{
    if (m_method==1)
    {
        if (code.compare("")==0) return 1;
        else if (code.compare("0")==0) return (t>0) ? 0 : -t;
        else if (code.compare("1")==0) return (t<0) ? 0 : t;
        else
        {
            double sup, inf, tcode = BinaryTree::getValueFromCode(code);
            computeBoundariesForBasisFunction(tcode,&inf,&sup,axis);
            if (t <= tcode && t >= inf) return ((t-inf)/(tcode-inf));
            else if (t >= tcode && t <= sup) return ((t-sup)/(tcode-sup));
            else return 0;
        }
    }
    else
    {
        if (code.compare("")==0) return 1;
        //else if (code.compare("0")==0) return (t>0) ? 0 : 0.5*t*(t - 1);
        //else if (code.compare("1")==0) return (t<0) ? 0 : 0.5*t*(t + 1);
        else if (code.compare("0")==0) return (t>0) ? 0 : -t*(t + 2);
        else if (code.compare("1")==0) return (t<0) ? 0 :  t*(2 - t);
        else
        {
            double sup, inf, tcode = BinaryTree::getValueFromCode(code);
            computeBoundariesForBasisFunction(tcode,&inf,&sup,axis);
            if (t <= sup && t >= inf) return (4/pow(sup-inf,2))*(sup-t)*(t-inf);
            else return 0;
        }
    }
}
/******************************************************************************/
