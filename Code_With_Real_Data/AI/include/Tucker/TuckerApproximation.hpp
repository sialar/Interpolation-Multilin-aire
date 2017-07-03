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

    private:
        static vector<string> keys;
        static int dimension;

        string core;
        string infoFileName;
        map<string,string> listOfFiles;
        map<string,vector<double>> listOfDomainBorders;
        map<string,vector<double>> listOfTuckerGridNodes;

        map<string,map<string,vector<vector<double>>>> finalOrthNormalizedEigVects;
        map<string,vector<double>> FinalTuckerCoeffs;
        map<string,vector<vector<int>>> listOfFinalCoefIndexes_arr;
        map<string,map<string,vector<vector<LagrangePolynomial>>>> listOfBasicFctsUsingLagrangeInterpolation;

    public:
        vector<string> listOfCrossSectionNames;

        ~TuckerApproximation() {};
        TuckerApproximation(string _core, vector<string> listOfCrossSection);

        static double computeKinf(map<string, double> listOfValues);
        static double find_nearest(vector<double> myList, double value);
        static double getInterpolation(vector<LagrangePolynomial> fct, double x);
        static vector<double> getInterpolationArr(vector<LagrangePolynomial> fct, MultiVariatePoint<double> x);
        static string check_string(string strs_begin, string strs_end, string NameFile);
        static vector<string> get_list_check_string(string strs_begin, string strs_end, string NameFile);
        static string convert_multiLines_oneLine(string lines);

        static vector<vector<double>> convert_str_dic_eigVects(string strs, string strs_split);
        static map<string,vector<double>> convert_str_dic_map(string strs, string strs_split);
        static vector<double> convert_str_dic_tucker_coef(string strs, string strs_split);
        static vector<vector<int>> convert_str_dic_coef_index(string strs, string strs_split);

        vector<vector<LagrangePolynomial>> getListOfInterpolationFcts(string Axis_k, string csName);
        double evaluate(MultiVariatePoint<double> point, string csName);
        map<string,vector<vector<double>>> getFinalOrthNormalizedEigVects(string NameFile);
        vector<double> getFinalTuckerCoeffs(string NameFile);
        vector<vector<int>> getListOfFinalCoefIndexes_arr(string NameFile);
        map<string,vector<vector<LagrangePolynomial>>> getListOfBasicFcts(string NameFile);

        void setListOfCrossSection(vector<string> listOfCrossSection);
        void setTuckerGridNodes();
        void setDomainBorders();
        void setListOfFiles();
        void setAll();
};

typedef std::shared_ptr<TuckerApproximation> TuckerApproximationPtr;
#endif
