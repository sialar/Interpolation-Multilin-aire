#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;
string Utils::projectPath = "/home/sialar/Stage/LaboJ_LLions/Code/";

void Utils::separateur()
{
    for (int i=0; i<206; i++) cout << "=";
    cout << endl;
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

double Utils::randomValue(double a, double b)
{
  return (rand()/(double)RAND_MAX) * (b-a) + a;
}

MultiVariatePoint<double> Utils::createRandomMultiVariatePoint(int d)
{
    MultiVariatePoint<double> point(d,0,0);
    for (int i=0; i<d; i++)
        point(i) = Utils::randomValue(0,1);
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
    ofstream file(projectPath + "data/dichotomy_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
        file << length << endl;
        vector<double> seq = createSequenceByDichotomy(length);
        for (int i=0; i<length; ++i)
              file << seq[i] << endl;
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::storeLejaSequenceInFile(int length)
{
    ofstream file(projectPath + "data/leja_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
        vector<double> lejaSeq = createLejaSequence(length);
        for (int i=0; i<length; ++i)
              file << lejaSeq[i] << endl;
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

vector<double> Utils::loadLejaSequenceFromFile(int length)
{
    ifstream file(projectPath + "data/leja_sequence.txt", ios::in);
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
        cerr << "Erreur à l'ouverture du fichier! " << endl;
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

int Utils::chooseDimensionD(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension d: ";
        cin >> dim;
    }
    return dim;
}

int Utils::chooseDimensionN(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension n: ";
        cin >> dim;
    }
    return dim;
}

int Utils::chooseNbTestPoints(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int nbTestPoints = -1;
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}

int Utils::chooseMaxIteration(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int maxIteration = -1;
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  Utils::separateur();
  return maxIteration;
}

int Utils::chooseMethod(int argc, char* argv[], int argNum)
{
    int method = -1;
    if (argc > argNum) method = stoi(argv[argNum]);
    while (method!=1 && method!=2)
    {
        cout << " - Choose the method of interpolation: " << endl;
        cout << "\t - 1: Using piecewise functions and middle points: " << endl;
        cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
        cin >> method;
    }
    return method;
}

bool Utils::saveResults(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  char store = 'x';
  while (store!='y' && store!='n')
  {
      cout << " - Store path and interpolation progression? (y/n) ";
      cin >> store;
  }
  return (store=='y');
}

bool Utils::saveError(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  char e = 'x';
  while (e!='y' && e!='n')
  {
      cout << " - Store interpolation error at the end of the algorithm? (y/n) ";
      cin >> e;
  }
  return (e=='y');
}

bool Utils::plotPath(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  char plot = 'x';
  while (plot!='y' && plot!='n')
  {
      cout << " - Save and plot interpolation points? (y/n) ";
      cin >> plot;
  }
  return (plot=='y');
}

int Utils::chooseFunction(int argc, char* argv[], int argNum)
{
  int f = 0;
  if (argc > argNum) f = stoi(argv[argNum]);
  while (f!=1 && f!=2 && f!=3 && f!=4)
  {
      cout << " - Choose the function to interpolate (1, 2, 3 or 4)";
      cin >> f;
  }
  return f;
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
