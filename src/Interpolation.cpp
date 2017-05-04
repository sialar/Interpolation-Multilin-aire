#include "../include/Interpolation.hpp"

Interpolation::~Interpolation()
{}

Interpolation::Interpolation(vector<int> sizes)
{
    m_d = sizes.size();
    m_points.resize(m_d);
    for (int i=0; i<m_d; i++)
        m_points[i].resize(sizes[i]);
}

void Interpolation::clear()
{
    m_path.clear();
    for (int i=0; i<m_d; i++)
        m_points[i].clear();
    m_points.clear();
    m_alphaMap.clear();
    m_curentNeighbours.clear();
}

/************************* Points *********************************************/
MultiVariatePoint<double> Interpolation::getPoint(MultiVariatePoint<int> nu)
{
    MultiVariatePoint<double> point(m_d,0.0);
    for (int i=0; i<m_d; i++)
        point(i) = m_points[i][nu(i)];
    return point;
}
void Interpolation::displayPoints()
{
    for (int i=0; i<m_d; i++)
    {
        cout << "Interpolation points in direction " << i << ": { ";;
        for (int j=0; j<int(m_points[i].size()); j++)
            cout << m_points[i][j] << " ";
        cout << "}" << endl;
    }
}
/******************************************************************************/

/************************* Path ***********************************************/
void Interpolation::savePathInFile()
{
    ofstream file("python/output_algo_AI.txt", ios::out | ios::trunc);
    if(file)
    {
        if (m_d==2)
        {
            cout << " - The path is saved in python/output_algo_AI.txt" << endl;
            file << m_points[0].size() << " " <<  m_points[1].size() << endl;
            file << m_path.size() << endl;
            for (MultiVariatePoint<int> nu : m_path)
                file << nu(0) << " " << nu(1) << endl;
        }
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}
void Interpolation::displayPath()
{
    cout << "Chemin = ";
    for (int i=0; i<int(m_path.size()); i++)
    {
        cout << "(";
        for (int k=0; k<m_d-1; k++)
            cout << m_path[i](k) << ",";
        cout << m_path[i](m_d-1) << ")";
        if (i != int(m_path.size()-1))
            cout << " --> ";
    }
    cout << endl;
}
bool Interpolation::indiceInPath(MultiVariatePoint<int>& index)
{
    for (MultiVariatePoint<int> nu : m_path)
        if (index==nu)
            return true;
    return false;
}
double Interpolation::tryWithCurentPath()
{
    vector<double> realValues, estimate;
    for (MultiVariatePoint<double> p : m_testPoints)
    {
        realValues.push_back(Utils::gNd(p));
        estimate.push_back(lagrangeInterpolation_ND_iterative(p));
    }
    return Utils::interpolationError(realValues,estimate);
}
int Interpolation::buildPathWithAIAlgo(int maxIteration, double start_time, double threshold)
{
    m_curentNeighbours.clear();
    m_path.clear();

    MultiVariatePoint<int> nu(m_d,0);
    m_path.push_back(nu);
    m_alphaMap[nu] = Utils::gNd(getPoint(nu));

    double val, max;
    MultiVariatePoint<int> argmax(m_d,0);
    int iteration = 1;

    while (iteration < maxIteration)
    {
        max = -numeric_limits<double>::max();
        updateCurentNeighbours(argmax);
        for (MultiVariatePoint<int> nu : m_curentNeighbours)
        {
            map<MultiVariatePoint<int>,double>::iterator it = m_alphaMap.find(nu);
            if (it != m_alphaMap.end())
                val = m_alphaMap[nu];
            else
                val =  computeLastAlphaNu(nu);


            if (val >= max)
            {
                max = val;
                argmax = nu;
            }
        }
        m_path.push_back(argmax);
        iteration++;

        // Test with curent path and evaluate the interpolation error on test points
        // If the error is lower than a threshold : stop AI
        if (iteration%(maxIteration/10)==0)
        {
            double run_time = omp_get_wtime() - start_time;
            double error = tryWithCurentPath();
            cout << "   - Interpolation error after " << iteration << " iterations: " << error;
            cout << " | Elapsed time : "  << run_time << endl;
            if (error < threshold)
            {
                cout << endl << "   - AI Algo stop after " << iteration << " iterations";
                cout << " | Elapsed time : "  << run_time << endl;
                return iteration;
            }
        }
    }
    return iteration;
}
double Interpolation::testPathBuilt(int maxIteration, double threshold)
// return the number of iteration
{
    double start_time = omp_get_wtime();
    int nbIterations = buildPathWithAIAlgo(maxIteration, start_time, threshold);
    double run_time = omp_get_wtime() - start_time;
    cout << endl << " - Time required to compute the path with " << nbIterations <<
            " iterations: " << run_time << "s" << endl;
    return run_time;
}
/******************************************************************************/

/************************* Alpha **********************************************/
void Interpolation::displayAlphaTab()
{
    map<MultiVariatePoint<int>,double>::iterator it;
    cout << "Values of alpha (" << m_alphaMap.size() << ") :" << endl;
    for (it=m_alphaMap.begin(); it!=m_alphaMap.end(); it++)
        cout << get<0>(*it) << ": " << get<1>(*it) << endl;
}
double Interpolation::computeLastAlphaNu(MultiVariatePoint<int>& nu)
{
    double lagrangeBasisFuncProd;
    double res = Utils::gNd(getPoint(nu));
    for (MultiVariatePoint<int> l : m_path)
    {
        lagrangeBasisFuncProd = 1;
        for (int p=0; p<m_d; p++)
            lagrangeBasisFuncProd *= lagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);
        res -= m_alphaMap[l] * lagrangeBasisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(nu,res));
    return res;
}
double Interpolation::testAlphaNuComputation(MultiVariatePoint<int>& nu)
{
    double start_time = omp_get_wtime();
    double val = computeLastAlphaNu(nu);
    double run_time = omp_get_wtime() - start_time;
    cout << " - Time required to compute alpha(" << nu << "): " << run_time << " s";
    cout << " Correct Value = " << val << endl;
    return val;
}
/******************************************************************************/

/************************* Neighbours *****************************************/
void Interpolation::displayCurentNeighbours()
{
    cout << "Curent neighbours (" << m_curentNeighbours.size() << ") = ";
    for (MultiVariatePoint<int> nu : m_curentNeighbours)
        cout << nu << " ";
    cout << endl;
}
bool Interpolation::isCorrectNeighbourToCurentPath(MultiVariatePoint<int>& nu)
{
    bool isNeighbour[m_d];
    for (int i=0; i<m_d; i++)
        isNeighbour[i] = true;

    MultiVariatePoint<int> index(nu);
    for (int i=0; i<m_d; i++)
        if (nu(i))
        {
            index(i)--;
            isNeighbour[i] = indiceInPath(index);
            index(i)++;
        }
    bool res = true;
    for (int i=0; i<m_d; i++)
        res = res && isNeighbour[i];
    return res;
}
void Interpolation::updateCurentNeighbours(MultiVariatePoint<int>& nu)
{
    m_curentNeighbours.remove(nu);
    MultiVariatePoint<int> nu1(nu);
    int max_indice[m_d];
    for (int i=0; i<m_d; i++)
        max_indice[i] = m_points[i].size()-1;
    for (int k=0; k<nu.getD(); k++)
    {
        if (nu1(k)<max_indice[k])
        {
            nu1(k)++;
            if (isCorrectNeighbourToCurentPath(nu1))
                m_curentNeighbours.push_back(nu1);
            nu1(k)--;
        }
    }
}
/******************************************************************************/


/*********************** Interpolation ****************************************/
void Interpolation::displayInterpolationPoints()
{
    cout << " - Dimension d = " << m_d << endl;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_points[i].size() << " points in direction " << i << " : { ";
        for (size_t j=0; j<m_points[i].size()-1; ++j)
            cout << m_points[i][j] << " ; ";
        cout << m_points[i][m_points[i].size()-1] << " }" << endl;
    }
}

double Interpolation::lagrangeBasisFunction_1D(int j, int k, double t, int axis)
{
    if (!k) return 1;
    double prod = 1;
    for (int i=0; i<k; ++i)
        if (i!=j)
            prod *= (t-m_points[axis][i]) / (m_points[axis][j]-m_points[axis][i]) ;
    return prod;
}
double Interpolation::lagrangeInterpolation_ND_iterative(MultiVariatePoint<double> x)
{
    double l_prod, sum = 0;
    for (MultiVariatePoint<int> nu : m_path)
    {
        l_prod = 1;
        for (int i=0; i<m_d; i++)
            l_prod *= lagrangeBasisFunction_1D(nu(i),nu(i),x(i),i);
        sum += l_prod * m_alphaMap[nu];
    }
    return sum;
}
double Interpolation::computeExecTimeOfOneApprox()
{
    MultiVariatePoint<double> x(m_d,0);
    for (int i=0; i<m_d; i++)
        x(i) = Utils::randomValue(-1,1);

    double start_time = omp_get_wtime();
    lagrangeInterpolation_ND_iterative(x);
    double run_time = omp_get_wtime() - start_time;
    return run_time;
}
/******************************************************************************/
