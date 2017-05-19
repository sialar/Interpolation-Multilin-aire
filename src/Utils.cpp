#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;

void Utils::separateur()
{
    cout << endl;
    for (int i=0; i<206; i++)
        cout << "*";
    cout << endl << endl;
}

void Utils::displayPoints(vector<double> points)
{
    cout << "{ ";
    for (size_t i=0; i<points.size()-1; ++i)
        cout << points[i] << " ; ";
    cout << points[points.size()-1] << " }" << endl;
}

void Utils::displayPoints(vector<MultiVariatePoint<double>*> points)
{
    cout << "{ ";
    for (size_t i=0; i<points.size()-1; ++i)
        cout << (*points[i]) << " ; ";
    cout << (*points[points.size()-1]) << " }" << endl;
}

double Utils::randomValue(double a, double b)
{
  return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

MultiVariatePoint<double>* Utils::createRandomMultiVariatePoint(int d)
{
    MultiVariatePoint<double>* point = new MultiVariatePoint<double>(d,0);
    for (int i=0; i<d; i++)
        (*point)(i) = Utils::randomValue(-1,1);
    return point;
}

double Utils::interpolationError(vector<double> realValue, vector<double> estimate)
{
    double e = 0;
    int n = min(int(realValue.size()),int(estimate.size()));
    for (int k=0; k<n; k++)
        e += abs(estimate[k] - realValue[k]);
    return e/n;
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

void Utils::storeResult(vector<MultiVariatePoint<double>*> x, vector<double> approx, vector<double> realValue)
{
    ofstream file("python/interpolation_result.txt", ios::out | ios::trunc);
    if(file)
    {
        for (int i=0; i<int(x.size()); i++)
        {
            for (int d=0; d<x[i]->getD(); d++)
                file << (*x[i])(d) << " ";
            cout << approx[i] << " " << realValue[i] << endl;
        }
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::storeFunction(vector<double> x, vector<double> y, vector<double> z)
{
  ofstream file("python/test.txt", ios::out | ios::trunc);
  if(file)
  {
      file << x.size() << endl;
      for (int i=0; i<int(x.size()); i++)
      {
          file << x[i] << " " << y[i] << " " << z[i] << endl;
      }
      file.close();
  }
  else
      cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::storeLejaSequenceInFile(vector<vector<double>> x)
{
    ofstream file("python/leja_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
      for (size_t i=0; i<x.size(); ++i)
      {
          file << x[i].size() << endl;
          for (size_t j=0; j<x[i].size(); ++j)
              file << x[i][j] << endl;
      }
      file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::storeLejaSequenceInFile(int length)
{
    ofstream file("lejaSequence", ios::out | ios::trunc);
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
    ifstream file("lejaSequence", ios::in);
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
        cerr << "Erreur à l'ouverture du fichier!" << endl;
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
double Utils::gNd(MultiVariatePoint<double>* x)
{
    double temp = 0;
    for (int i=0; i<x->getD(); i++)
    temp += pow((*x)(i),2);
    return sin(sqrt(temp));

    //double temp = 0;
    //for (int i=0; i<x->getD(); i++)
    //temp += pow((*x)(i),2);
    //return temp;

    //double temp = 0;
    //for (int i=0; i<x->getD(); i++)
    //  temp += sqrt(1 - pow((*x)(i),2));
    //return temp;
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
