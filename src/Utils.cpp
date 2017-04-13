#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;


void Utils::storeLejaSequenceInFile(vector<double> seq)
{
    ofstream file("leja/leja_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
        file << seq.size() << endl;
        for (double dx : seq)
        {
            for (double dy : seq)
            {
                file << dx << " " << dy << endl;
            }
        }
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
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
    m_1dGrid = createChebychevSequence(10000);
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
    storeLejaSequenceInFile(points);
    return points;
}

bool Utils::isTooCloseToOneLejiPoint(double y, vector<double> seq, double threshold)
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
        if (!isTooCloseToOneLejiPoint(m_1dGrid[k],seq,1e-10))
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

double Utils::randomValue(double a, double b)
{
    return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

double Utils::squareError(vector<double> realValue, vector<double> estimate)
{
    double e = 0;
    if (realValue.size()!=estimate.size())
        cerr << "Erreur: les deux tableaux ne sont pas de la même taille!" << endl;
    int n = min(int(realValue.size()),int(estimate.size()));
    for (int k=0; k<n; k++)
        e += pow(estimate[k] - realValue[k],2);
    return e/n;
}

void Utils::separateur()
{
    cout << endl;
    for (int i=0; i<100; i++)
        cout << "*";
    cout << endl << endl;
}
