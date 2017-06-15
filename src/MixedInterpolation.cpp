#include "../include/MixedInterpolation.hpp"

MixedInterpolation::MixedInterpolation(int d, int n, int nIter, MultiVariatePoint<int> methods, Function f) :
    Interpolation(d,n,nIter,f)
{
    m_lejaSequence = Utils::loadLejaSequenceFromFile(m_maxIteration);
    m_trees.resize(m_d);
    m_methods = methods;
    for (int i=0; i<m_d; i++)
      m_trees[i] = make_shared<BinaryTree>();
}
void MixedInterpolation::clearAllTrees()
{
    for (int i=0; i<int(m_trees.size()); i++)
        m_trees[i]->clearTree();
}

/************************* Data Points ****************************************/
MultiVariatePoint<double> MixedInterpolation::getPoint(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<double> point(m_d,0,0.0);
    for (int i=0; i<m_d; i++)
    {
        if (m_methods(i))
            point(i) = BinaryTree::getValueFromCode((*nu)(i));
        else
            point(i) = m_lejaSequence[stoi((*nu)(i))];
    }
    return point;
}
void MixedInterpolation::addInterpolationPoint(MultiVariatePoint<double>p)
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
            if (m_methods(i)) m_trees[i]->addNode(p(i));
        }
    }
}
void MixedInterpolation::computeBoundariesForBasisFunction(double t, double* inf, double* sup, int axis)
{
    /*
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
    */
    // Or using trees
    m_trees[axis]->searchNode(t,inf,sup,false);
}
/******************************************************************************/


/************************* AI algo ********************************************/
bool customAlphaLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return Utils::norm(nu->getAlpha(),2) < Utils::norm(mu->getAlpha(),2);
}
bool customAgeLess(MultiVariatePointPtr<string> nu, MultiVariatePointPtr<string> mu)
{
    return nu->getWaitingTime() < mu->getWaitingTime();
}
MultiVariatePointPtr<string> MixedInterpolation::maxElement(int iteration)
{
    if (iteration%4)
        return *max_element(m_curentNeighbours.begin(),m_curentNeighbours.end(),customAlphaLess);
    else
    {
        MultiVariatePointPtr<string> mu = *max_element(m_curentNeighbours.begin(), \
                                            m_curentNeighbours.end(),customAgeLess);
        for (MultiVariatePointPtr<string> nu : m_curentNeighbours)
            if (nu->alphaIsNull() && nu->getWaitingTime()==mu->getWaitingTime())
                return nu;
        return mu;
    }
}
MultiVariatePointPtr<string> MixedInterpolation::getFirstMultivariatePoint()
{
    MultiVariatePointPtr<string> nu = make_shared<MultiVariatePoint<string>>(m_d,m_n,"");
    for (int i=0; i<m_d; i++)
        if (!m_methods(i))
            (*nu)(i) = "0";
    return nu;
}
void MixedInterpolation::updateCurentNeighbours(MultiVariatePointPtr<string> nu)
{
    m_curentNeighbours.remove(nu);
    for (MultiVariatePointPtr<string> mu : m_curentNeighbours)
        mu->incrWaitingTime();

    vector<string> childrenCodes;
    for (int i=0; i<m_d; i++)
    {
        if (m_methods(i))
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
        else
        {
            MultiVariatePointPtr<string> mu = make_shared<MultiVariatePoint<string>>(*nu);
            (*mu)(i) = to_string(stoi((*mu)(i))+1);
            if (isCorrectNeighbourToCurentPath(mu))
                m_curentNeighbours.push_back(mu);
        }
    }
}

bool MixedInterpolation::isCorrectNeighbourToCurentPath(MultiVariatePointPtr<string> nu)
{
    MultiVariatePoint<string> mu(*nu);
    string temp;
    for (int i=0; i<m_d; i++)
    {
        if (m_methods(i))
        {
            if ((*nu)(i).compare("")!=0)
            {
                temp = mu(i);
                mu(i) = BinaryTree::getParentCode(temp);
                if (!indiceInPath(mu)) return false;
                mu(i) = temp;
            }
        }
        else
        {
            if ((*nu)(i).compare("0")!=0)
            {
                mu(i) = to_string(stoi(mu(i))-1);
                if (!indiceInPath(mu)) return false;
                mu(i) = to_string(stoi(mu(i))+1);
            }
        }
    }
    return true;
}
bool MixedInterpolation::indiceInPath(MultiVariatePoint<string> index)
{
    for (MultiVariatePointPtr<string> nu : m_path)
    if (Utils::equals(index,*nu))
        return true;
    return false;
}
/******************************************************************************/

/************************* Display functions **********************************/
void MixedInterpolation::savePathInFile(string fileName)
{
  ofstream file(fileName, ios::out | ios::trunc);
  if(file)
  {
    if (m_d==2 || m_d==3)
    {
        for (int i=0; i<m_d; i++)
            file << m_interpolationPoints[i].size() << " ";
        file << endl;
        for (int i=0; i<m_d; i++)
            file << m_methods(i) << " ";
        file << endl << m_path.size() << endl;
        MultiVariatePoint<double> x(m_d,0,0.0);
        for (MultiVariatePointPtr<string> nu : m_path)
        {
            x = getPoint(nu);
            for (int i=0; i<m_d; i++)
                file << x(i) << " ";
            for (int i=0; i<m_d; i++)
                file << (*nu)(i) << " " ;
            file << endl;
        }
        for (int i=0; i<int(m_path.size())-1; i++)
            file << Utils::vector2str(m_path[i]->getAlpha()) << " ; ";
        file << Utils::vector2str(m_path[m_path.size()-1]->getAlpha()) << endl;
    }
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
void MixedInterpolation::saveInterpolationBasisFunctions()
{
  ofstream file(Utils::projectPath + "data/basis_functions.txt", ios::out | ios::trunc);
  if(file)
  {
    if (m_d==1)
    {
      vector<double> x;
      int nbPoints = int(m_testPoints.size());
      for (int i=0; i<nbPoints; i++)
          x.push_back(m_testPoints[i](0));
      sort(x.begin(), x.end());
      file << m_path.size() << " " << nbPoints << " " << m_methods(0) << endl;
      MultiVariatePoint<double> p;
      for (int j=0; j<nbPoints; j++)
      {
        file << x[j];
        for (MultiVariatePointPtr<string> nu : m_path)
            file << " " << basisFunction_1D((*nu)(0),x[j],0);
        p = MultiVariatePoint<double>::toMonoVariatePoint(x[j]);
        file << " " <<  Utils::vector2str(func(p));
        file << endl;
      }
      for (int i=0; i<int(m_path.size())-1; i++)
          file << Utils::vector2str(m_path[i]->getAlpha()) << " ; ";
      file << Utils::vector2str(m_path[m_path.size()-1]->getAlpha()) << endl;
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
vector<double> MixedInterpolation::tryWithDifferentMethods(MultiVariatePoint<int> methods, double threshold)
{
    clearAll();
    clearAllTrees();
    cout << "   - Interpolation using methods " << methods;
    vector<double> errors;
    setMethods(methods);
    testPathBuilt(threshold, m_maxIteration<21);
    vector<vector<double>> realValues, estimate;
    for (MultiVariatePoint<double> p : m_testPoints)
    {
        realValues.push_back(func(p));
        estimate.push_back(interpolation(p,m_path.size()));
    }
    errors.push_back(Utils::relativeInterpolationError(realValues,estimate));
    errors.push_back(Utils::mseInterpolationError(realValues,estimate));
    return errors;
}
bool customLess(vector<double> e1, vector<double> e2)
{
    return Utils::norm(e1,2) < Utils::norm(e2,2);
}
bool relativeErrorLess(vector<double> e1, vector<double> e2)
{
  return e1[0] < e2[0];
}
bool mseErrorLess(vector<double> e1, vector<double> e2)
{
    return e1[1] < e2[1];
}
MultiVariatePoint<int> MixedInterpolation::tryAllCases(double threshold)
{
    MultiVariatePoint<int> methods(m_d, 0, 0);
    vector<vector<double>> errors(3);
    for (int i=0; i<3; i++) errors[i].resize(2,0.0);
    vector<double> temp(2,0.0);
    for (int i=0; i<m_d; i++)
    {
        Utils::separateur();
        cout << " - Comparing methods in direction " << i << ": " << endl;
        for (int j=0; j<nbMethods; j++)
        {
            methods(i) = j;
            map<MultiVariatePoint<int>,vector<double>>::iterator it = m_methods_errors.find(methods);
            if (it == m_methods_errors.end())
            {
                temp = tryWithDifferentMethods(methods, threshold);
                errors[j][0] = temp[0];
                errors[j][1] = temp[1];
                m_methods_errors.insert(pair<MultiVariatePoint<int>,vector<double>>(methods,errors[j]));
            }
            else errors[j] = get<1>(*it);
            cout << "   ---> Relative Interpolation error when using " << methods << " = " << errors[j][0] << endl;
            cout << "   ---> MSE Interpolation error when using " << methods << " = " << errors[j][1] << endl;
        }
        methods(i) = distance(errors.begin(),min_element(errors.begin(), errors.end(),customLess));
        cout << endl << " ---> Chosen method in direction " << i << ": " <<  methods(i) << endl;
    }
    Utils::separateur();
    cout << " - The optimal choice of methods is " << methods << endl;
    cout << " - Relative Interpolation error (pcm) = " << (*min_element(errors.begin(), errors.end(), relativeErrorLess))[0] << endl;
    cout << " - MSE Interpolation error (pcm) = " << (*min_element(errors.begin(), errors.end(), mseErrorLess))[1] << endl;
    if (Utils::norm(*min_element(errors.begin(), errors.end(), customLess),2) > 1)
    {
        cout << " - The interpolation error is higher than 1 pcm, ";
        cout << "you should increase the maximum number of iterations" << endl;
    }
    return methods;
}
double MixedInterpolation::basisFunction_1D(string code, double t, int axis)
{
    if (m_methods(axis)==0)
    {
        if (code.compare("0")==0) return 1;
        int k = stoi(code);
        double prod = 1;
        for (int i=0; i<k; ++i)
            prod *= (t-m_lejaSequence[i]) / (m_lejaSequence[k]-m_lejaSequence[i]) ;
        return prod;
    }
    else if (m_methods(axis)==1)
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
