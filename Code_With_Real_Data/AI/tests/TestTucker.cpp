#include <python3.5/Python.h>
#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    string core = "MOX/";
    vector<string> listOfCrossSectionNamesWithArgs;
    listOfCrossSectionNamesWithArgs.push_back("macro_totale0");
    listOfCrossSectionNamesWithArgs.push_back("macro_totale1");
    listOfCrossSectionNamesWithArgs.push_back("macro_absorption0");
    listOfCrossSectionNamesWithArgs.push_back("macro_absorption1");
    listOfCrossSectionNamesWithArgs.push_back("macro_scattering000");
    listOfCrossSectionNamesWithArgs.push_back("macro_scattering001");
    listOfCrossSectionNamesWithArgs.push_back("macro_scattering010");
    listOfCrossSectionNamesWithArgs.push_back("macro_scattering011");
    listOfCrossSectionNamesWithArgs.push_back("macro_nu*fission0");
    listOfCrossSectionNamesWithArgs.push_back("macro_nu*fission1");
    listOfCrossSectionNamesWithArgs.push_back("macro_fission0");
    listOfCrossSectionNamesWithArgs.push_back("macro_fission1");

    map<string, string> listOfFiles;
    listOfFiles.insert(pair<string,string>("macro_totale0",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_totale0.txt"));
    listOfFiles.insert(pair<string,string>("macro_totale1",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_totale1.txt"));
    listOfFiles.insert(pair<string,string>("macro_absorption0",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_absorption0.txt"));
    listOfFiles.insert(pair<string,string>("macro_absorption1",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_absorption1.txt"));
    listOfFiles.insert(pair<string,string>("macro_scattering000",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_scattering000.txt"));
    listOfFiles.insert(pair<string,string>("macro_scattering001",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_scattering001.txt"));
    listOfFiles.insert(pair<string,string>("macro_scattering010",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_scattering010.txt"));
    listOfFiles.insert(pair<string,string>("macro_scattering011",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_scattering011.txt"));
    listOfFiles.insert(pair<string,string>("macro_nu*fission0",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_nu_fission0.txt"));
    listOfFiles.insert(pair<string,string>("macro_nu*fission1",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_nu_fission1.txt"));
    listOfFiles.insert(pair<string,string>("macro_fission0",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_fission0.txt"));
    listOfFiles.insert(pair<string,string>("macro_fission1",Utils::projectPath+"AI/data/"+core+"TuckerInfor_macro_fission1.txt"));


    vector<int> listOfSubdivisionDirection(5);
    vector<vector<double>> listOfValuesForSubdivision(5);
    vector<vector<int>> listOfNumberOfPointsForSubdivision(5);

    vector<double> listOfValuesForSubdivision_0;
    listOfValuesForSubdivision_0.push_back(150);
    listOfValuesForSubdivision_0.push_back(10000);
    vector<int> listOfNumberOfPointsForSubdivision_0;
    listOfNumberOfPointsForSubdivision_0.push_back(3);
    listOfNumberOfPointsForSubdivision_0.push_back(9);
    listOfSubdivisionDirection[0] = 1;
    listOfValuesForSubdivision[0] = listOfValuesForSubdivision_0;
    listOfNumberOfPointsForSubdivision[0] = listOfNumberOfPointsForSubdivision_0;

    vector<double> listOfValuesForSubdivision_2;
    listOfValuesForSubdivision_2.push_back(0.71692097187);
    vector<int> listOfNumberOfPointsForSubdivision_2;
    listOfNumberOfPointsForSubdivision_2.push_back(5);
    listOfSubdivisionDirection[2] = 1;
    listOfValuesForSubdivision[2] = listOfValuesForSubdivision_2;
    listOfNumberOfPointsForSubdivision[2] = listOfNumberOfPointsForSubdivision_2;

    string NameFile = Utils::projectPath + "AI/data/" + core + "GeneralInfor_MOX.txt";
    string strs_begin = "listOfDomainBorders = {";
    string strs_end = "}";

    string lines = TuckerApproximation::check_string(strs_begin, strs_end, NameFile);
    string converted_line = TuckerApproximation::convert_multiLines_oneLine(lines);

    string strs_split = " = ";
    vector<vector<double>> listOfDomainBorders = TuckerApproximation::convert_str_dic(converted_line, strs_split);
    strs_begin = "listOfTuckerGridNodes = {";
    strs_end = "}";

    lines = TuckerApproximation::check_string(strs_begin, strs_end, NameFile);

    converted_line = TuckerApproximation::convert_multiLines_oneLine(lines);
    string str_old = "array(";
    string str_new = "";
    converted_line = TuckerApproximation::replace_str(converted_line, str_old, str_new);

    str_old = ")";
    str_new = "";
    converted_line = TuckerApproximation::replace_str(converted_line, str_old, str_new);

    strs_split = " = ";
    vector<vector<double>> listOfTuckerGridNodes_dic = TuckerApproximation::convert_str_dic( converted_line, strs_split);

    vector<vector<double>> listOfTuckerGridNodes = listOfTuckerGridNodes_dic;

    map<string,vector<vector<vector<double>>>> finalOrthNormalizedEigVects;
    map<string,vector<double>> FinalTuckerCoeffs;
    map<string,vector<vector<int>>> listOfFinalCoefIndexes_arr;
    map<string,vector<vector<vector<LagrangePolynomial>>>> listOfBasicFctsUsingLagrangeInterpolation;

    for (string csName : listOfCrossSectionNamesWithArgs)
    {
        NameFile = listOfFiles[csName];
        finalOrthNormalizedEigVects.insert(pair<string,vector<vector<vector<double>>>>(csName, \
                                            TuckerApproximation::getFinalOrthNormalizedEigVects(NameFile)));
        FinalTuckerCoeffs.insert(pair<string,vector<double>>(csName, TuckerApproximation::getFinalTuckerCoeffs(NameFile)));
        listOfFinalCoefIndexes_arr.insert(pair<string,vector<vector<int>>>(csName, TuckerApproximation::getListOfFinalCoefIndexes_arr(NameFile)));
        listOfBasicFctsUsingLagrangeInterpolation.insert(pair<string,vector<vector<vector<LagrangePolynomial>>>>( \
                                csName, TuckerApproximation::getListOfBasicFcts(listOfDomainBorders, listOfTuckerGridNodes, \
                                finalOrthNormalizedEigVects[csName])));
    }
    vector<double> point(5,0.1);
    double value_T;
    vector<double> values;
    vector<string> listOfCrossSectionNames = listOfCrossSectionNamesWithArgs;
    for (string csName : listOfCrossSectionNames)
    {
        value_T = TuckerApproximation::evaluate(listOfBasicFctsUsingLagrangeInterpolation[csName], \
                  listOfFinalCoefIndexes_arr[csName], FinalTuckerCoeffs[csName], point);
        values.push_back(value_T);
        cout << value_T << " ";
    }
    cout << endl;
    return 0;
}
