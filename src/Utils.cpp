#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;

void Utils::storeResult1D(vector<double> x, vector<double> y, vector<double> real_y)
{
    ofstream file("python/interpolation_result_1D.txt", ios::out | ios::trunc);
    if(file)
    {
        file << x.size() << endl;
        for (int i=0; i<int(x.size()); i++)
            file << x[i] << " " << y[i] << " " << real_y[i] << endl;
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}


void Utils::storeResult2D(vector<double> x, vector<double> y, vector<double> z, vector<double> real)
{
    ofstream file("python/interpolation_result_2D.txt", ios::out | ios::trunc);
    if(file)
    {
        file << x.size() << " " << y.size() << endl;
        for (int i=0; i<int(x.size()); i++)
            for (int j=0; j<int(y.size()); j++)
                file << x[i] << " " << y[j] << " " << z[int(y.size())*i+j] <<
                         " " << real[int(y.size())*i+j] << endl;
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::store2DLejaSequenceInFile(vector<double> x, vector<double> y)
{
    ofstream file("python/leja_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
        file << x.size()<< endl;
        for (int i=0; i<int(x.size()); i++)
            file << x[i] << endl;
        file << y.size()<< endl;
        for (int i=0; i<int(y.size()); i++)
            file << y[i] << endl;
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Utils::store3DLejaSequenceInFile(vector<double> x, vector<double> y, vector<double> z)
{
    ofstream file("python/leja_sequence.txt", ios::out | ios::trunc);
    if(file)
    {
        file << x.size()<< endl;
        for (int i=0; i<int(x.size()); i++)
            file << x[i] << endl;
        file << y.size()<< endl;
        for (int i=0; i<int(y.size()); i++)
            file << y[i] << endl;
        file << z.size()<< endl;
        for (int i=0; i<int(z.size()); i++)
            file << z[i] << endl;
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
    for (int i=0; i<206; i++)
        cout << "*";
    cout << endl << endl;
}

double Utils::g1d(double y)
{
    return y*y*y + pow((y+1)*sin(y),3) + y*exp(2*y) + 1;
}

double Utils::g2d(double x, double y)
{
    //return sin(sqrt(x*x + y*y));
    return x*x*exp(y) + pow(y*sin(2*x*x),3);
    //return x*x + y*y*x + 2*y + 1;
}

double Utils::gNd(vector<double> x)
{
    double temp = 0;
    for (int i=0; i<int(x.size()); i++)
        temp += pow(x[i],2);
    return sin(sqrt(temp));
}

void Utils::displayPoints(vector<vector<double>> v, int d)
{
    for (int k=0; k<d; k++)
    {
        cout << "         + Suivant la direction " << k << ": ";
        for (int i=0; i<int(v[k].size()); ++i)
            cout << v[k][i] << " ";
        cout << endl;
    }
}

vector<double> Utils::displayGRealValues(vector<double> vx, vector<double> vy, int d, bool debug)
{
    vector<double> realValues;
    double val = 0;
    if (d==1)
    {
        cout << "   - Calcul direct (evaluation de g en " << vx.size() << " points):" << endl;
        for (int i=0; i<int(vx.size()); i++)
        {
            val = Utils::g1d(vx[i]);
            realValues.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }
    else if (d==2)
    {
        cout << "   - Calcul direct (evaluation de g en " << vx.size()*vy.size() << " points):" << endl;
        for (int i=0; i<int(vx.size()); i++)
        {
            for (int j=0; j<int(vy.size()); j++)
            {
                val = Utils::g2d(vx[i],vy[j]);
                realValues.push_back(val);
                if (debug) cout << val << " ";
            }
            if (debug) cout << endl;
        }
    }
    return realValues;
}

void Utils::displayApproximation(vector<double> approx, int nx, int ny, int d, bool debug)
{
    if (debug)
    {
        for (int i=0; i<nx; i++)
        {
            if (d==1)
                cout << approx[i] << " ";
            else if (d==2)
            {
                for (int j=0; j<ny; j++)
                    cout << approx[i*nx+j] << " ";
                cout << endl;
            }
        }
        cout << endl;
    }
}
