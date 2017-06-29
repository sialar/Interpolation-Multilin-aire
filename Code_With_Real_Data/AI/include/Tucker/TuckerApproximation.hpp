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

using namespace std;

class TuckerApproximation
{
    public:
        static vector<string> keys;
        static int dimension;

        static double getInterpolation(vector<LagrangePolynomial> fct, double x);
        static vector<double> getInterpolationArr(vector<LagrangePolynomial> fct, vector<double> x);

        static double computeKinf(map<string, double> listOfValues);
        static double find_nearest(vector<double> myList, double value);
        static vector<vector<LagrangePolynomial>> getListOfInterpolationFcts( string Axis_k, \
                                      map<string,vector<double>> listOfDomainBorders, \
                                      map<string,vector<double>> listOfTuckerGridNodes, \
                                      map<string,vector<vector<double>>> finalOrthNormalizedEigVects);

        static double evaluate(map<string,vector<vector<LagrangePolynomial>>> listOfBasisFcts, \
                               vector<vector<int>> listOfFinalCoefIndexes_arr, \
                               vector<double> FinalTuckerCoeffs, vector<double> point);

        static string check_string(string strs_begin, string strs_end, string NameFile);
        static vector<string> get_list_check_string(string strs_begin, string strs_end, string NameFile);
        static string convert_multiLines_oneLine(string lines);
        static string replace_str(string strs, string str_old, string str_new);
        static vector<vector<double>> convert_str_dic_eigVects(string strs, string strs_split);
        static map<string,vector<double>> convert_str_dic_map(string strs, string strs_split);
        static vector<double> convert_str_dic_tucker_coef(string strs, string strs_split);
        static vector<vector<int>> convert_str_dic_coef_index(string strs, string strs_split);
        static map<string,vector<vector<double>>> getFinalOrthNormalizedEigVects(string NameFile);
        static vector<double> getFinalTuckerCoeffs(string NameFile);
        static vector<vector<int>> getListOfFinalCoefIndexes_arr(string NameFile);
        static map<string,vector<vector<LagrangePolynomial>>> getListOfBasicFcts( \
                                    map<string,vector<double>> listOfDomainBorders, \
                                    map<string,vector<double>>  listOfTuckerGridNodes, \
                                    map<string,vector<vector<double>>> finalOrthNormalizedEigVects);
};

#endif
