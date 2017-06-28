#ifndef TUCKERAPPROXIMATION
#define TUCKERAPPROXIMATION

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <memory>
#include <cmath>
#include <chrono>

#include "LagrangePolynomial.hpp"

using namespace std;

class TuckerApproximation
{

    public:
        static int dimension;

        static double getInterpolation(vector<LagrangePolynomial> fct, double x);
        static vector<double> getInterpolationArr(vector<LagrangePolynomial> fct, vector<double> x);

        static double computeKinf(map<string, double> listOfValues);
        static double find_nearest(vector<double> myList, double value);
        static vector<vector<LagrangePolynomial>> getListOfInterpolationFcts( int Axis_k, \
                                      vector<vector<double>> listOfDomainBorders, \
                                      vector<vector<double>> listOfTuckerGridNodes, \
                                      vector<vector<vector<double>>> finalOrthNormalizedEigVects);

        static double evaluate(vector<vector<vector<LagrangePolynomial>>> istOfBasisFcts, \
                               vector<vector<int>> listOfFinalCoefIndexes_arr, \
                               vector<double> FinalTuckerCoeffs, vector<double> point);

        static string check_string(string strs_begin, string strs_end, string NameFile);
        static vector<string> get_list_check_string(string strs_begin, string strs_end, string NameFile);
        static string convert_multiLines_oneLine(string lines);
        static string replace_str(string strs, string str_old, string str_new);
        static vector<vector<double>> convert_str_dic(string strs, string strs_split);
        static vector<vector<int>> convert_str_dic_int(string strs, string strs_split);
        static vector<vector<vector<double>>> getFinalOrthNormalizedEigVects(string NameFile);
        static vector<double> getFinalTuckerCoeffs(string NameFile);
        static vector<vector<int>> getListOfFinalCoefIndexes_arr(string NameFile);
        static vector<vector<vector<LagrangePolynomial>>> getListOfBasicFcts( \
                                    vector<vector<double>> listOfDomainBorders, \
                                    vector<vector<double>>  listOfTuckerGridNodes, \
                                    vector<vector<vector<double>>> finalOrthNormalizedEigVects);
};

#endif
