#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int chooseDimension(int argc, char* argv[])
{
    if (argc > 1) return stoi(argv[1]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension : ";
        cin >> dim;
    }
    return dim;
}

int chooseNbTestPoints(int argc, char* argv[])
{
  if (argc > 2) return stoi(argv[2]);
  int nbTestPoints = -1;
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}

int chooseMaxIteration(int argc, char* argv[])
{
  if (argc > 3) return stoi(argv[3]);
  int maxIteration = -1;
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  Utils::separateur();
  return maxIteration;
}

void saveAllErrorsInFile(vector<MultiVariatePoint<double>> errors, vector<int> iterations)
{
  ofstream file("data/interpolation_error.txt", ios::out | ios::trunc);
  if(file)
  {
      file << errors.size() << endl;
      int size = errors.size();
      for (int i=0; i<size; i++)
          file << iterations[i] << " " << errors[i](0) << " " << errors[i](1) \
               << " " << errors[i](2) << endl;
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    Utils::separateur();
    int dim = chooseDimension(argc,argv);
    int nbTestPoints = chooseNbTestPoints(argc,argv);
    int maxIteration = chooseMaxIteration(argc,argv);
    // Initialisation of test points
    vector<MultiVariatePoint<double>> testPoints;
    vector<double> realValues, estimate;
    testPoints.resize(nbTestPoints);
    for (int j=0; j<nbTestPoints; j++)
        testPoints[j] = Utils::createRandomMultiVariatePoint(dim);

    // Method 1
    LagrangeInterpolationPtr interp_0(new LagrangeInterpolation(dim,maxIteration,Utils::g));
    interp_0->setTestPoints(testPoints);
    interp_0->setSaveError(true);
    // Path creation
    double threshold = 1e-10;
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_0->testPathBuilt(threshold, maxIteration<21);
    // Computing real values, and approximation of function g at test points
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(interp_0->func(p));
        estimate.push_back(interp_0->interpolation_ND(p,interp_0->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    realValues.clear();
    estimate.clear();
    Utils::separateur();

    // Method 2
    PiecewiseInterpolationPtr interp_1(new PiecewiseInterpolation(dim,maxIteration,1,Utils::g));
    interp_1->setTestPoints(testPoints);
    interp_1->setSaveError(true);
    // Path creation
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_1->testPathBuilt(threshold, maxIteration<21);
    // Computing real values, and approximation of function g at test points
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(interp_1->func(p));
        estimate.push_back(interp_1->interpolation_ND(p,interp_1->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    realValues.clear();
    estimate.clear();
    Utils::separateur();


    // Method 3
    PiecewiseInterpolationPtr interp_2(new PiecewiseInterpolation(dim,maxIteration,2,Utils::g));
    interp_2->setTestPoints(testPoints);
    interp_2->setSaveError(true);
    // Path creation
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_2->testPathBuilt(threshold, maxIteration<21);
    // Computing real values, and approximation of function g at test points
    for (MultiVariatePoint<double> p : testPoints)
    {
        realValues.push_back(interp_2->func(p));
        estimate.push_back(interp_2->interpolation_ND(p,interp_2->path().size()));
    }
    cout << " - Interpolation error = " << Utils::interpolationError(realValues,estimate) << endl;
    Utils::separateur();


    vector<MultiVariatePoint<double>> allErrors;
    map<int,double>::iterator it;
    vector<int> iterations;
    double epsilon = 1e-10;
    for (it=interp_0->errors().begin();it!=interp_0->errors().end(); it++)
        iterations.push_back(get<0>(*it));
    for (int i=0; i<int(iterations.size()); i++)
    {
        MultiVariatePoint<double> p(3,0);
        p(0) = interp_0->errors()[iterations[i]];
        p(1) = interp_1->errors()[iterations[i]];
        p(2) = interp_2->errors()[iterations[i]];
        if (p(0)>epsilon && p(1)>epsilon && p(2)>epsilon)
            allErrors.push_back(p);
    }
    saveAllErrorsInFile(allErrors, iterations);

    return 0;
}
