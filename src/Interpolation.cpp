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
void Interpolation::buildPathWithAIAlgo(int k, bool parallel)
{
    m_curentNeighbours.clear();
    m_path.clear();

    MultiVariatePoint<int> nu(m_d,0);
    m_path.push_back(nu);
    m_alphaMap[nu] = Utils::gNd(getPoint(nu));

    double val, max;
    MultiVariatePoint<int> argmax(m_d,0);
    int iteration = 1;

    while (iteration < k)
    {
        max = -numeric_limits<double>::max();
        updateCurentNeighbours(argmax);
        for (MultiVariatePoint<int> nu : m_curentNeighbours)
        {
            map<MultiVariatePoint<int>,double>::iterator it = m_alphaMap.find(nu);
            if (it != m_alphaMap.end())
                val = m_alphaMap[nu];
            else if (parallel)
                val =  computeLastAlphaNuPar(nu);
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
    }
}
double Interpolation::testPathBuilt(int nbIteration, bool parallel)
{
    std::clock_t chrono_work_start = std::clock();
    timeval chrono_elpased_start; gettimeofday(&chrono_elpased_start, NULL);

    buildPathWithAIAlgo(nbIteration, parallel);

    timeval chrono_elpased_stop; gettimeofday(&chrono_elpased_stop, NULL);
    std::clock_t chrono_work_stop  = std::clock();
    double elapsedTimeSeconds = double(chrono_elpased_stop.tv_sec - chrono_elpased_start.tv_sec)
                            + double(chrono_elpased_stop.tv_usec - chrono_elpased_start.tv_usec)*1e-6  ;
    double totalWorkSeconds = double(chrono_work_stop - chrono_work_start) / CLOCKS_PER_SEC;
    cout << " - Time required to compute the path with " << nbIteration <<
            " iterations: " << elapsedTimeSeconds << " s";
    cout << " (total work = " << totalWorkSeconds << " s)" << endl;
    return elapsedTimeSeconds;
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
double Interpolation::computeLastAlphaNuPar(MultiVariatePoint<int>& nu)
{
    double res = Utils::gNd(getPoint(nu));
    #pragma omp parallel for reduction(+:res) schedule(guided)
    for (int i=0; i<int(m_path.size()); i++)
    {
        double lagrangeBasisFuncProd = 1;
        MultiVariatePoint<int> l = m_path[i];
        for (int p=0; p<m_d; p++)
            lagrangeBasisFuncProd *= lagrangeBasisFunction_1D(l(p),l(p),m_points[p][nu(p)],p);

        res -= m_alphaMap[l] * lagrangeBasisFuncProd;
    }
    m_alphaMap.insert(pair<MultiVariatePoint<int>,double>(nu,res));
    return res;
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
    std::clock_t chrono_work_start = std::clock();
    timeval chrono_elpased_start; gettimeofday(&chrono_elpased_start, NULL);

    double val = computeLastAlphaNu(nu);

    timeval chrono_elpased_stop; gettimeofday(&chrono_elpased_stop, NULL);
    std::clock_t chrono_work_stop  = std::clock();
    double elapsedTimeSeconds = double(chrono_elpased_stop.tv_sec - chrono_elpased_start.tv_sec)
                            + double(chrono_elpased_stop.tv_usec - chrono_elpased_start.tv_usec)*1e-6  ;
    double totalWorkSeconds = double(chrono_work_stop - chrono_work_start) / CLOCKS_PER_SEC;
    cout << " - Time required to compute alpha(" << nu << "): " << elapsedTimeSeconds << " s";
    cout << " (total work = " << totalWorkSeconds << "s)";
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

    float temps;
    clock_t t1, t2;
    t1 = clock();
    lagrangeInterpolation_ND_iterative(x);
    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    return temps;
}
/******************************************************************************/
