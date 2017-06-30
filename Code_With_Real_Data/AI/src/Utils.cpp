#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;
string Utils::projectPath = "/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/";

void Utils::separateur()
{
    for (int i=0; i<206; i++) cout << "=";
    cout << endl;
}

void Utils::displayValues(vector<double> values)
{
    cout << "[ ";
    for (size_t i=0; i<values.size()-1; ++i)
        cout << values[i] << " ; ";
    cout << values[values.size()-1] << " ]" << endl;
}

void Utils::displayPoints(vector<MultiVariatePoint<double>> points)
{
    cout << "{ ";
    for (size_t i=0; i<points.size()-1; ++i)
        cout << points[i] << " ; ";
    cout << points[points.size()-1] << " }" << endl;
}

void Utils::displayPoints(vector<vector<double>> points)
{
    int n = points[0].size();
    cout << "{ ";
    for (size_t i=0; i<points.size()-1; ++i)
    {
        cout << "(";
        for (int j=0; j<n-1; ++j)
            cout << points[i][j] << ",";
        cout << points[i][n-1] << ")";
        cout << " ; ";
    }
    cout << "(";
    for (size_t i=0; i<points[points.size()-1].size()-1; ++i)
        cout << points[points.size()-1][i] << ",";
    cout << points[points.size()-1][n-1] << ")";
    cout << " }" << endl;
}

double Utils::adaptCoordsToFunctionDomain(double a, double b, double x)
{
    return (2*x)/(b-a) + 1 - (2*b)/(b-a);
}

double Utils::convertToFunctionDomain(double a, double b, double x)
{
    return x*(b-a)/2 + (a+b)/2;
}

MultiVariatePoint<double> Utils::getCoordsFromString(string s)
{
    MultiVariatePoint<double> coords(5,0,0.0);
    s.pop_back();
    s.erase(0,1);
    stringstream ss(s);
    string subs;
    int axis = 0;
    while (getline(ss, subs, ','))
    {
        coords(axis) = double(stod(subs));
        axis++;
    }
    return coords;
}

double Utils::randomValue(double a, double b)
{
  return (rand()/(double)RAND_MAX) * (b-a) + a;
}

MultiVariatePoint<double> Utils::createRandomMultiVariatePoint(int d)
{
    MultiVariatePoint<double> point(d,0,0);
    for (int i=0; i<d; i++)
        point(i) = Utils::randomValue(-1,1);
    return point;
}

double Utils::relativeInterpolationError(vector<vector<double>> realValue, vector<vector<double>> estimate)
{
    vector<double> e, f;
    int n = min(int(realValue.size()),int(estimate.size()));
    for (int k=0; k<n; k++)
    {
        e.push_back(norm(diff(estimate[k],realValue[k]),0));
        f.push_back(norm(realValue[k],0));
    }
    return pow(10,5) * (*max_element(e.begin(), e.end()) / *max_element(f.begin(), f.end()));
}

double Utils::mseInterpolationError(vector<vector<double>> realValue, vector<vector<double>> estimate)
{
    double e = 0.0, f = 0.0;
    int n = min(int(realValue.size()),int(estimate.size()));
    for (int k=0; k<n; k++)
    {
        e += pow(norm(diff(estimate[k],realValue[k]),0),2);
        f += pow(norm(realValue[k],0),2);
    }
    return pow(10,5) * (sqrt(e) / sqrt(f));
}

bool normLess(double x, double y)
{
    return abs(x) < abs(y);
}
double Utils::norm(vector<double> x, int p)
{
    if (!p) return abs(*max_element(x.begin(), x.end(), normLess));
    else
    {
        double res = 0.0;
        for (int i=0; i<int(x.size()); i++)
            res += pow(x[i],2);
        return sqrt(res);
    }
}

vector<double> Utils::diff(vector<double> x,vector<double> y)
{
    vector<double> res;
    for (int i=0; i<int(x.size()); i++)
        res.push_back(x[i]-y[i]);
    return res;
}

bool Utils::isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold)
{
    for (int i=0; i<int(seq.size()); i++)
        if (abs(y-seq[i]) <= threshold)
            return true;
    return false;
}

double Utils::computeNewLejaPointFromSequence(vector<double> seq)
{
    int argmax = -1;
    double prod = 1;
    double max = numeric_limits<double>::min();
    for (int k=0; k<int(m_1dGrid.size()); k++)
    {
        prod = 1;
        if (!isTooCloseToOneLejaPoint(m_1dGrid[k],seq,1e-10))
        {
            for (int i=0; i<int(seq.size()); i++)
                prod *= abs(m_1dGrid[k] - seq[i]);
            if (prod > max)
            {
                max = prod;
                argmax = k;
            }
        }
    }
    return m_1dGrid[argmax];
}

void Utils::binaryDecomposition(int number, vector<double>& binary_decomp)
{
  binary_decomp.clear();
  int temp = number;
  while (temp > 0)
  {
    binary_decomp.push_back(temp % 2);
    temp = temp / 2;
  }
}

void Utils::storeDichotomySequenceInFile(int length)
{
    ofstream file(projectPath + "AI/data/dichotomy_sequence.dat", ios::out | ios::trunc);
    if(file)
    {
        file << length << endl;
        vector<double> seq = createSequenceByDichotomy(length);
        for (int i=0; i<length; ++i)
              file << seq[i] << endl;
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}

void Utils::storeLejaSequenceInFile(int length)
{
    ofstream file(projectPath + "AI/data/leja_sequence.dat", ios::out | ios::trunc);
    if(file)
    {
        vector<double> lejaSeq = createLejaSequence(length);
        for (int i=0; i<length; ++i)
              file << lejaSeq[i] << endl;
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}

vector<double> Utils::loadLejaSequenceFromFile(int length)
{
    ifstream file(projectPath + "AI/data/leja_sequence.dat", ios::in);
    vector<double> lejaSeq;
    string line;
    int line_index = 0;
    if(file)
    {
        while (getline(file,line) && line_index<length)
        {
            lejaSeq.push_back(stof(line));
            line_index++;
        }
        file.close();
    }
    else
        cerr << "Error while opening the file! " << endl;
    return lejaSeq;
}

vector<double> Utils::createSequenceByDichotomy(int length)
{
    vector<double> sequence;
    if (length > 0) sequence.push_back(0);
    if (length > 1) sequence.push_back(-1);
    if (length > 2) sequence.push_back(1);
    int e = 0, k = 3;
    double denom;
    while (k<length)
    {
        e++;
        denom = pow(2,e);
        for (int j=denom-1; j>=0 && int(sequence.size())<length; j--)
            if (j%2) sequence.push_back(-j/denom);
        for (int j=0; j<denom && int(sequence.size())<length; j++)
            if (j%2) sequence.push_back(j/denom);
        k += denom;
    }
    return sequence;
}

string Utils::vector2str(vector<double> x)
{
    string s = "(" + to_string(x[0]);
    for (int i=1; i<int(x.size()); i++)
        s = s + ", " + to_string(x[i]);
    s = s + ") ";
    return s;
}

vector<double> Utils::createUniformSequence(int nbPoints)
{
    double sum = 0;
    vector<double> points;
    vector<double> binary_decomp;
    int stop = (nbPoints%2) ? (nbPoints-1)/2 : (nbPoints-1)/2+1;
    points.resize(nbPoints);

    points[0] =  1;
    points[1] = -1;
    points[2] =  0;
    for (int k=1; k<stop; ++k)
    {
        sum = 0;
        binaryDecomposition(k,binary_decomp);
        for (size_t j=0; j<binary_decomp.size(); ++j)
            sum += binary_decomp[j] / pow(2,j);

        points[2*k+1] = 0.5 * sum;
        points[2*k+2] = - points[2*k+1];
    }
    return points;
}

vector<double> Utils::createChebychevSequence(int nbPoints)
{
    vector<double> points;
    points.resize(nbPoints);
    for (int i=0; i<nbPoints; i++)
        points[i] = cos(i*M_PI/(nbPoints-1));
    return points;
}

vector<double> Utils::createLejaSequence(int nbPoints)
{
    m_1dGrid = createChebychevSequence(10001);
    vector<double> points;
    points.push_back(1);
    double newPoint;
    for (int i=1; i<nbPoints; i++)
    {
        while (int(points.size()) < nbPoints)
        {
            newPoint = computeNewLejaPointFromSequence(points);
            points.push_back(newPoint);
        }
    }
    return points;
}

bool Utils::equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2)
{
    if (nu1.getD() != nu2.getD()) return false;
    for (int i=0; i<nu1.getD(); i++)
        if (nu1(i).compare(nu2(i))!=0)
            return false;
    return true;
}

bool Utils::displayResults()
{
  char display = 'x';
  while (display!='y' && display!='n')
  {
      cout << " - Display path and interpolation points: (y/n) " ;
      cin >> display;
  }
  return (display=='y');
}
