#ifndef TUCKERAPPROXIMATION
#define TUCKERAPPROXIMATION

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>


#include "LagrangePolynomial.hpp"
#include "../MultiVariatePoint.hpp"
#include "../Utils.hpp"

using namespace std;

class TuckerApproximation
{

    public:
        static vector<string> keys;
        static int dimension;

        string core;
        string infoFileName;
        vector<string> listOfCrossSectionNames;
        map<string,string> listOfFiles;
        map<string,vector<double>> listOfDomainBorders;
        map<string,vector<double>> listOfTuckerGridNodes;

        map<string,map<string,vector<vector<double>>>> finalOrthNormalizedEigVects;
        map<string,vector<double>> FinalTuckerCoeffs;
        map<string,vector<vector<int>>> listOfFinalCoefIndexes_arr;
        map<string,map<string,vector<vector<LagrangePolynomial>>>> listOfBasicFctsUsingLagrangeInterpolation;

        ~TuckerApproximation() {};
        TuckerApproximation(string _core, vector<string> listOfCrossSection);
        TuckerApproximation(string _core, string crossSection);

        static double computeKinf(map<string, double> listOfValues);
        static double find_nearest(vector<double> myList, double value);

        void setListOfCrossSection(vector<string> listOfCrossSection);
        void setTuckerGridNodes();
        void setDomainBorders();
        void setListOfFiles();
        void setAll();

        map<string,vector<vector<double>>> getFinalOrthNormalizedEigVects(string NameFile);
        vector<double> getFinalTuckerCoeffs(string NameFile);
        vector<vector<int>> getListOfFinalCoefIndexes_arr(string NameFile);
        map<string,vector<vector<LagrangePolynomial>>> getListOfBasicFcts(string csName);

        vector<vector<LagrangePolynomial>> getListOfInterpolationFcts(string Axis_k, string csName);
        double evaluate(MultiVariatePoint<double> point, string csName);
        double getInterpolation(vector<LagrangePolynomial> fct, double x);


};

typedef std::shared_ptr<TuckerApproximation> TuckerApproximationPtr;
#endif
