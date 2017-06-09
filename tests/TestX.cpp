#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include "../include/Utils.hpp"
#include "../include/BinaryTree.hpp"
#include "../include/MultiVariatePoint.hpp"
#include "../include/Interpolation.hpp"
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/LagrangeInterpolation.hpp"

using namespace std;

MultiVariatePoint<double> extractInfoFromLine(string line)
{
    MultiVariatePoint<double> data;
    //TODO
    return data;
}

vector<MultiVariatePoint<double>> parseFile(string fileName)
{
    vector<MultiVariatePoint<double>> data;
    ifstream file(fileName, ios::in);
    if(file)
    {
        string line;
        while(getline(file, line))
        {
            cout << line << endl;
            MultiVariatePoint<double> pointInfo = extractInfoFromLine(line);
            data.push_back(pointInfo);
        }
        file.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    return data;
}

double evaluateCrossSection(vector<MultiVariatePoint<double>> dataPoints,
                            vector<double> dataValues, int method,
                            MultiVariatePoint<double>& x, int iteration, Function f)
{
    int dim = x.getD();
    if (method)
    {
        PiecewiseInterpolationPtr interp(new PiecewiseInterpolation(dim,1,pow(10,dim),method,f));
        interp->testPathBuilt(1e-5, false, 0);
        return interp->interpolation(x,interp->path().size());
    }
    else
    {
        LagrangeInterpolationPtr interp(new LagrangeInterpolation(dim,1,pow(10,dim),f));
        interp->testPathBuilt(1e-5, false, 0);
        return interp->interpolation(x,interp->path().size());
    }
}


int main( int argc, char* argv[] )
{
    vector<MultiVariatePoint<double>> physicalParamaters;
    vector<double> crossSection;

    physicalParamaters = parseFile("../data/data.txt");

    return 0;
}
