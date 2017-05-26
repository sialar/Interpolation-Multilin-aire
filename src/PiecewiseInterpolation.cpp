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
int PiecewiseInterpolation::buildPathWithAIAlgo(auto start_time, double threshold, bool debug)
{
    m_curentNeighbours.clear();
    m_path.clear();
    MultiVariatePointPtr<string> argmax = make_shared<MultiVariatePoint<string>>(m_d,"");
    m_curentNeighbours.push_back(argmax);
    int iteration = 0;
    double val;

    while (!m_curentNeighbours.empty() && iteration < m_maxIteration)
    {
        if (debug)
        {
            Utils::separateur();
            displayPath();
            displayCurentNeighbours();
        }
        for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
        {
            if (nu->alphaAlreadyComputed()) val = nu->getAlpha();
            else val =  computeLastAlphaNu(nu);
            if (debug) cout << *nu << " " << val << " | ";
        }
        if (debug) cout << endl;
        argmax = *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),(iteration%100) ? alphaLess : ageLess);
        m_path.push_back(argmax);
        addInterpolationPoint(getPoint(argmax));
        updateCurentNeighbours(argmax);
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        /*
        if ((m_maxIteration>10) && iteration%(m_maxIteration/10)==0)
        {
            auto end_time = chrono::steady_clock::now();
            std::chrono::duration<double> run_time = end_time - start_time;
            double error = tryWithCurentPath();
            cout << "   - Interpolation error after " << iteration << " iterations: " << error;
            cout << " | Elapsed time : "  << run_time.count() << endl;
            if (error < threshold)
            {
                cout << endl << "   - AI Algo stop after " << iteration << " iterations";
                cout << " | Elapsed time : "  << run_time.count() << endl;
                return iteration;
            }
        }
        */

    }
    return iteration;
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
double PiecewiseInterpolation::computeLastAlphaNu(MultiVariatePointPtr<string> nu)
{
    double basisFuncProd;
    double res = Utils::gNd(getPoint(nu));
    for (MultiVariatePointPtr<string> l : m_path)
    {
        basisFuncProd = 1;
        for (int p=0; p<m_d; p++)
        {
            if (m_method == 1)
                basisFuncProd *= piecewiseFunction_1D((*l)(p),getPoint(nu)(p),p);
            else if (m_method == 2)
                basisFuncProd *= quadraticFunction_1D((*l)(p),getPoint(nu)(p),p);
        }
        res -= l->getAlpha() * basisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePointPtr<string>,double>(nu,res));
    nu->setAlpha(res);
    return res;
}
/******************************************************************************/


/*********************** Test functions ***************************************/
void PiecewiseInterpolation::testPathBuilt(double threshold, bool debug)
{
  auto start_time = chrono::steady_clock::now();
  int nbIterations = buildPathWithAIAlgo(start_time, threshold, debug);
  auto end_time = chrono::steady_clock::now();
  std::chrono::duration<double> run_time = end_time - start_time;
  cout << endl << " - Time required to compute the path with " << nbIterations <<
  " iterations: " << run_time.count() << "s" << endl;
}
double PiecewiseInterpolation::tryWithCurentPath()
{
  vector<double> realValues, estimate;
  for (MultiVariatePoint<double> p : m_testPoints)
  {
    realValues.push_back(Utils::gNd(p));
    estimate.push_back(interpolation_ND(p,m_path.size()));
  }
  return Utils::interpolationError(realValues,estimate);
}
/******************************************************************************/


/************************* Display functions **********************************/
void PiecewiseInterpolation::displayPath()
{
    // Format: (nu:points[nu]) --> () ...
    cout << "Chemin =";
    int n = m_path.size();
    for (int i=0; i<n; i++)
    {
        if (i>0) cout << "\t";
        cout << " " << i << " :";
        cout << " [" << *m_path[i] << ":" << getPoint(m_path[i]) << "]";
        cout << " --> alpha" << *m_path[i] << " = " << m_path[i]->getAlpha() << endl;
    }
    cout << endl;
}
void PiecewiseInterpolation::displayAlphaTab()
{
    map<MultiVariatePointPtr<string>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << "(" << *(get<0>(*it)) << ":" << getPoint((get<0>(*it))) << ") : " << get<1>(*it) <<  endl;
}
void PiecewiseInterpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
        cout << "(" << (*nu) << ":" << getPoint(nu) << ":" << nu << ") [" << nu->getWaitingTime() << "] | ";
    cout << endl << endl;
}

void PiecewiseInterpolation::displayTrees()
{
    for (int i=0; i<m_d; i++)
    {
        cout << " - Binary Tree in direction " << i << " :" << endl;
        m_trees[i]->displayBinaryTree();
        cout << endl;
    }
}

void PiecewiseInterpolation::savePathInFile()
{
  ofstream file("data/path.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2)
    {
      //cout << " - The path is saved in data/path.txt" << endl;
      file << m_interpolationPoints[0].size() << " " <<  m_interpolationPoints[1].size() << endl;
      file << m_path.size() << endl;
      MultiVariatePoint<double> x(m_d,0.0);
      for (MultiVariatePointPtr<string> nu : m_path)
      {
        x = getPoint(nu);
        file << x(0) << " " << x(1) << " ";
        file << (*nu)(0) << " " << (*nu)(1) << endl;
      }
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
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
        {
          if (m_method == 1)
          file << " " << piecewiseFunction_1D((*nu)(0),x[j],0);
          else if (m_method == 2)
          file << " " << quadraticFunction_1D((*nu)(0),x[j],0);
        }
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

void PiecewiseInterpolation::storeInterpolationProgression()
{
  ofstream file("data/interpolation_progression.txt", ios::out | ios::trunc);
  if(file)
  {
      if (m_d==1)
      {
          vector<double> x;
          int nbPoints = int(m_testPoints.size());
          for (int i=0; i<nbPoints; i++)
              x.push_back(m_testPoints[i](0));
          sort(x.begin(), x.end());
          MultiVariatePoint<double> p;
          vector<double> tempPath;
          for (double t : x)
          {
              for (int i=0; i<int(m_path.size()); i++)
              {
                  p = MultiVariatePoint<double>::toMonoVariatePoint(t);
                  file << interpolation_ND(p, i+1) << " ";
              }
              file << endl;
          }
      }
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}

/******************************************************************************/


/*********************** Interpolation ****************************************/
double PiecewiseInterpolation::piecewiseFunction_1D(string code, double t, int axis)
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
double PiecewiseInterpolation::quadraticFunction_1D(string code, double t, int axis)
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
double PiecewiseInterpolation::interpolation_ND(MultiVariatePoint<double>& x, int end)
{
  double l_prod, sum = 0;
  for (int k=0; k<end; k++)
  {
      l_prod = 1;
      for (int i=0; i<m_d; i++)
      {
          if (m_method == 1)
              l_prod *= piecewiseFunction_1D((*m_path[k])(i),x(i),i);
          else if (m_method == 2)
              l_prod *= quadraticFunction_1D((*m_path[k])(i),x(i),i);
      }
      sum += l_prod * m_path[k]->getAlpha();
  }
  return sum;
}
/******************************************************************************/
