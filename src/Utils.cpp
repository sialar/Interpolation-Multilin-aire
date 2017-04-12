#include "../include/Utils.hpp"


vector<double> Utils::createChebychevSequence(int nbPoints)
{
    // http://www.uvt.rnu.tn/resources-uvt/cours/analyse_num/chap4/node8.html
    vector<double> points;
    points.resize(nbPoints);
    for (int i=0; i<nbPoints; i++)
        points[i] = cos((2*i+1)*M_PI/(2*nbPoints));
    return points;
}

vector<double> Utils::createLejaSequence(int nbPoints)
{
    double sum = 0;
    vector<double> points;
    vector<double> binary_decomp;
    points.resize(nbPoints);

    points[0] =  1;
    for (int k=1; k<nbPoints; ++k)
    {
        sum = 0;
        binaryDecomposition(k,binary_decomp);
        for (size_t j=0; j<binary_decomp.size(); ++j)
            sum += binary_decomp[j] / pow(2,j);

        points[k] = cos(M_PI*sum);
    }
    return points;
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
        cout << "Erreur: les deux tableaux ne sont pas de la même taille!" << endl;
    int n = min(int(realValue.size()),int(estimate.size()));
    for (int k=0; k<n; k++)
        e += pow(estimate[k] - realValue[k],2);
    return e/n;
}
