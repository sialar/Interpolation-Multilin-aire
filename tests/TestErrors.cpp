#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/LagrangeInterpolation.hpp"
#include "../include/PiecewiseInterpolation.hpp"
#include "../include/Utils.hpp"

using namespace std;


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
    int dimD = Utils::chooseDimensionD(argc,argv,1);
    int dimN = Utils::chooseDimensionN(argc,argv,2);
    int nbTestPoints = Utils::chooseNbTestPoints(argc,argv,3);
    int maxIteration = Utils::chooseMaxIteration(argc,argv,4);

    double threshold = 1e-20;
    vector<vector<double>> realValues, estimate;

    // Method 1
    LagrangeInterpolationPtr interp_0(new LagrangeInterpolation(dimD,dimN,maxIteration,Functions::f));
    interp_0->setRandomTestPoints(nbTestPoints);
    interp_0->setSaveError(true);
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_0->testPathBuilt(threshold, maxIteration<21);
    for (MultiVariatePoint<double> p : interp_0->testPoints())
    {
        realValues.push_back(interp_0->func(p));
        estimate.push_back(interp_0->interpolation(p,interp_0->path().size()));
    }
    cout << " - Relative Interpolation error (pcm) = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error (pcm) = " << Utils::mseInterpolationError(realValues,estimate) << endl;
    realValues.clear();
    estimate.clear();
    Utils::separateur();

    // Method 2
    PiecewiseInterpolationPtr interp_1(new PiecewiseInterpolation(dimD,dimN,maxIteration,1,Functions::f));
    interp_1->setRandomTestPoints(nbTestPoints);
    interp_1->setSaveError(true);
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_1->testPathBuilt(threshold, maxIteration<21);
    for (MultiVariatePoint<double> p : interp_1->testPoints())
    {
        realValues.push_back(interp_1->func(p));
        estimate.push_back(interp_1->interpolation(p,interp_1->path().size()));
    }
    cout << " - Relative Interpolation error (pcm) = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error (pcm) = " << Utils::mseInterpolationError(realValues,estimate) << endl;    realValues.clear();
    estimate.clear();
    Utils::separateur();


    // Method 3
    PiecewiseInterpolationPtr interp_2(new PiecewiseInterpolation(dimD,dimN,maxIteration,2,Functions::f));
    interp_2->setRandomTestPoints(nbTestPoints);
    interp_2->setSaveError(true);
    cout << " - The maximum number of iterations in AI algo: " << maxIteration << endl;
    interp_2->testPathBuilt(threshold, maxIteration<21);
    for (MultiVariatePoint<double> p : interp_2->testPoints())
    {
        realValues.push_back(interp_2->func(p));
        estimate.push_back(interp_2->interpolation(p,interp_2->path().size()));
    }
    cout << " - Relative Interpolation error (pcm) = " << Utils::relativeInterpolationError(realValues,estimate) << endl;
    cout << " - MSE Interpolation error (pcm) = " << Utils::mseInterpolationError(realValues,estimate) << endl;
    Utils::separateur();

    vector<MultiVariatePoint<double>> allErrors;
    map<int,double>::iterator it;
    vector<int> iterations;
    for (it=interp_0->errors().begin();it!=interp_0->errors().end(); it++)
        iterations.push_back(get<0>(*it));

    for (int i=0; i<int(iterations.size()); i++)
    {
        MultiVariatePoint<double> p(3,0,0);
        p(0) = interp_0->errors()[iterations[i]];
        p(1) = interp_1->errors()[iterations[i]];
        p(2) = interp_2->errors()[iterations[i]];
        if (p(0)>threshold && p(1)>threshold && p(2)>threshold)
            allErrors.push_back(p);
    }
    saveAllErrorsInFile(allErrors, iterations);
    return 0;
}
