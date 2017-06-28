#include "../../include/Tucker/TuckerApproximation.hpp"

int TuckerApproximation::dimension = 5;

double min_elt(vector<double> v)
{
    if (v.empty()) return -11;
    double m = v[0];
    for (double x : v)
        if (x<m)
            m = x;
    return m;
}

double max_elt(vector<double> v)
{
    if (v.empty()) return -11;
    double m = v[0];
    for (double x : v)
        if (x>m)
            m = x;
    return m;
}


double TuckerApproximation::getInterpolation(vector<LagrangePolynomial> fct, double x)
{
    int n = fct.size();
    double min = min_elt(fct[0].axisValues());
    double max = max_elt(fct[n-1].axisValues());

    cout << n << " " << min << " "  << x << " "   << max << " "   << endl;
    if (n == 1)  return fct[0].getInterpolation(x);
    if (x < min) return fct[0].getInterpolation(x);
    if (x > max) return fct[n-1].getInterpolation(x);

    if (n > 1 && min <= x && x <= max)
    {
        int j = 0;
        while (j <= n-1)
        {
            min = min_elt(fct[j].axisValues());
            max = max_elt(fct[j].axisValues());
            if (min <= x && x <= max)
            {
                return fct[j].getInterpolation(x);
                break;
            }
            else j++;
        }
    }
    return -12;

}

vector<double> TuckerApproximation::getInterpolationArr(vector<LagrangePolynomial> fct, vector<double> x)
{
    vector<double> listOfValues;
    for (double xi : x) listOfValues.push_back(getInterpolation(fct,xi));
    return listOfValues;
}

double TuckerApproximation::computeKinf(map<string, double> listOfValues)
{
    double nufi0 = listOfValues["macro_nu*fission0"];
    double nufi1 = listOfValues["macro_nu*fission1"];
    double tot0 = listOfValues["macro_totale0"];
    double tot1 = listOfValues["macro_totale1"];
    double scatt000 = listOfValues["macro_scattering000"];
    double scatt001 = listOfValues["macro_scattering001"];
    double scatt010 = listOfValues["macro_scattering010"];
    double scatt011 = listOfValues["macro_scattering011"];
    return ( nufi0 * (tot1 - scatt011) + nufi1*scatt010) / ( (tot0 - scatt000 ) * (tot1 - scatt011) - scatt001*scatt010);
}

double TuckerApproximation::find_nearest(vector<double> myList, double value)
{
    double min_dist = abs(myList[0]-value);
    for (double x : myList)
        if (abs(x-value)<min_dist)
            min_dist = abs(x-value);
    return min_dist;
}

vector<double> getWhereEquals(vector<double> vec, double val)
{
    vector<double> res;
    for (int i=0; i<int(vec.size()); i++)
        if (vec[i]==val) res.push_back(i);
    return res;
}

vector<double> vectorSequence(vector<double> vec, int i1, int i2)
{
    vector<double> res;
    for (int i=i1; i<i2; i++)
        res.push_back(vec[i]);
    return res;
}

vector<vector<LagrangePolynomial>> TuckerApproximation::getListOfInterpolationFcts( int Axis_k, \
                              vector<vector<double>> listOfDomainBorders, \
                              vector<vector<double>> listOfTuckerGridNodes, \
                              vector<vector<vector<double>>> finalOrthNormalizedEigVects)
{
    vector<vector<LagrangePolynomial>> listOfBasicFctsUsingLagrangeInterpolation_Axis_k;
    int nbExtremes =  listOfDomainBorders[Axis_k].size();
    int nbOfFcts_Axis_k = finalOrthNormalizedEigVects[Axis_k].size();
    for (int i=0; i<nbOfFcts_Axis_k; i++)
    {
        vector<LagrangePolynomial> polLagrange_ki;
        if (nbExtremes > 2)
        {
            for (int j=0; j<(nbExtremes - 1); j++)
            {
                double value1 = find_nearest(listOfTuckerGridNodes[Axis_k], listOfDomainBorders[Axis_k][j]);
                double value2 = find_nearest(listOfTuckerGridNodes[Axis_k], listOfDomainBorders[Axis_k][j+1]);

                vector<double> j1Arr = getWhereEquals(listOfTuckerGridNodes[Axis_k],value1);
                vector<double> j2Arr = getWhereEquals(listOfTuckerGridNodes[Axis_k],value2);

                int j1 = j1Arr.back();
                int j2 = j2Arr.front();
                j2 = j2 + 1;

                LagrangePolynomial interp_k_couplage(vectorSequence(listOfTuckerGridNodes[Axis_k],j1,j2),\
                                          vectorSequence(finalOrthNormalizedEigVects[Axis_k][i],j1,j2));

                polLagrange_ki.push_back(interp_k_couplage);
            }

        }
        else if  (nbExtremes == 2)
        {
            LagrangePolynomial interp(listOfTuckerGridNodes[Axis_k],finalOrthNormalizedEigVects[Axis_k][i]);
            polLagrange_ki.push_back(interp);
        }
        listOfBasicFctsUsingLagrangeInterpolation_Axis_k.push_back(polLagrange_ki);
    }
    return listOfBasicFctsUsingLagrangeInterpolation_Axis_k;
}

double TuckerApproximation::evaluate(vector<vector<vector<LagrangePolynomial>>> listOfBasisFcts, \
                                     vector<vector<int>> listOfFinalCoefIndexes_arr, \
                                     vector<double> FinalTuckerCoeffs, vector<double> point)
{
    double approxValue = 0.0;
    for (int i=0; i<int(listOfFinalCoefIndexes_arr.size()); i++)
    {
        double d = 1.0;
        vector<LagrangePolynomial> fct;
        for (int Axis_k=0; Axis_k<dimension; Axis_k++)
        {
            fct = listOfBasisFcts[Axis_k][listOfFinalCoefIndexes_arr[i][Axis_k]];
            d = d*getInterpolation(fct, point[Axis_k]);
        }
        approxValue = approxValue +  FinalTuckerCoeffs[i]*d;
    }
    return approxValue;
}

bool isSubString(string str, string substr)
{
    size_t found = str.find(substr);
    return (found!=string::npos);
}

string TuckerApproximation::check_string(string strs_begin, string strs_end, string NameFile)
{
    cout << NameFile << endl;
    ifstream file(NameFile, ios::in);
    if(file)
    {
        bool found = false;
        string found_strs = "";
        string line;
        while (getline(file, line))
        {
            if (isSubString(line,strs_begin) && isSubString(line,strs_end))
            {
                file.close();
                return line;
            }
            else if (isSubString(line,strs_begin) && !isSubString(line,strs_end))
            {
                found_strs = found_strs + line;
                while (!isSubString(line,strs_end))
                {
                    getline(file, line);
                    found_strs = found_strs + line;
                }
                found = true;
                file.close();
                return found_strs;
            }
        }
        if (!found) cout << "The translation cannot be found!" << endl;
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
    return "";
}


vector<string> TuckerApproximation::get_list_check_string(string strs_begin, string strs_end, string NameFile)
{
    vector<string> list_found_strs;
    ifstream file(NameFile, ios::in);
    if(file)
    {
        string line;
        while (getline(file, line))
        {
            string found_strs = "";
            if (isSubString(line,strs_begin) && !isSubString(line,strs_end))
            {
                found_strs = found_strs + line;
                while (isSubString(line,strs_end))
                {
                    getline(file, line);
                    found_strs = found_strs + line;
                }
                list_found_strs.push_back(found_strs);
            }
        }
        file.close();
    }
    else cerr << "Error while opening the file!" << endl;
    return list_found_strs;
}

string TuckerApproximation::convert_multiLines_oneLine(string lines)
{
    string converted_line = "";
    for (int i=0; i<int(lines.length()); i++)
        if (lines[i] != '\0') converted_line += lines[i];
        else converted_line += ' ';
    return converted_line;
}

string replace(string* str, string oldstr, string newstr)
{
    size_t found = str->find(oldstr);
    str->replace(found,oldstr.length(),newstr);
    return *str;
}

string TuckerApproximation::replace_str(string strs, string str_old, string str_new)
{
    if (isSubString(strs,str_old))
    {
        string new_strs = replace(&strs, str_old, str_new);
        return new_strs;
    }
    else
    {
      cout << strs << endl;
      cout << str_old << endl;
      cout << str_new << endl;
        cerr << "Warning : strs can not be replaced!" << endl;
        exit(1);
    }
    return "";
}

vector<vector<double>> TuckerApproximation::convert_str_dic(string strs, string strs_split)
{
    vector<vector<double>> dic;
    if (isSubString(strs,strs_split))
    {
        size_t found = strs.find(strs_split);
        string strs_dic =  strs.substr(found,strs.length());
        //dic = ast.literal_eval(strs_dic) //TODO
        return dic;
    }
    else cerr << "Warning : strs_split is not in strs!" << endl;
    return dic;
}

vector<vector<int>> TuckerApproximation::convert_str_dic_int(string strs, string strs_split)
{
    vector<vector<int>> dic;
    if (isSubString(strs,strs_split))
    {
        size_t found = strs.find(strs_split);
        string strs_dic =  strs.substr(found,strs.length());
        //dic = ast.literal_eval(strs_dic) //TODO
        return dic;
    }
    else cerr << "Warning : strs_split is not in strs!" << endl;
    return dic;
}

vector<vector<vector<double>>> TuckerApproximation::getFinalOrthNormalizedEigVects(string NameFile)
{
    string strs_begin = "Orthonormalized eigenvectors";
    string strs_end = "]]";
    vector<string> list_check_string = get_list_check_string(strs_begin, strs_end, NameFile);
    vector<vector<vector<double>>> finalOrthNormalizedEigVects;
    for (int i=0; i<dimension; i++)
    {
        string lines = list_check_string[i];
        string converted_line = convert_multiLines_oneLine(lines);
        string str_old = "[ ";
        string str_new = "[";
        if (isSubString(converted_line,str_old))
            converted_line = replace_str(converted_line, str_old, str_new);
        str_old = "]] ";
        str_new = "]]";
        converted_line = replace_str(converted_line, str_old, str_new);
        str_old = " ";
        str_new = ", ";
        converted_line = replace_str(converted_line, str_old, str_new);
        string strs_split = "self.finalOrthNormalizedEigVects_Axis_k:, ";
        finalOrthNormalizedEigVects.push_back(convert_str_dic(converted_line, strs_split));
    }
    return finalOrthNormalizedEigVects;
}

vector<double> TuckerApproximation::getFinalTuckerCoeffs(string NameFile)
{
    string strs_begin = "Coefficient values";
    string strs_end = "]";
    string lines = check_string(strs_begin, strs_end, NameFile);
    string converted_line = convert_multiLines_oneLine(lines);
    string str_old = "[ ";
    string str_new = "[";
    if (isSubString(converted_line,str_old))
        converted_line = replace_str(converted_line, str_old, str_new);

    str_old = "]: ";
    str_new = "]";
    if (isSubString(converted_line,str_old))
        converted_line = replace_str(converted_line, str_old, str_new);

    str_old = " ";
    str_new = ", ";
    converted_line = replace_str(converted_line, str_old, str_new);

    string strs_split = "self.FinalTuckerCoeffs, =";
    vector<double> FinalTuckerCoeffs;// = convert_str_dic(converted_line, strs_split); //TODO
    return FinalTuckerCoeffs;
}


vector<vector<int>> TuckerApproximation::getListOfFinalCoefIndexes_arr(string NameFile)
{
    string strs_begin = "self.listOfFinalCoefIndexes_arr";
    string strs_end = "]]";
    string lines = check_string(strs_begin, strs_end, NameFile);
    string converted_line = convert_multiLines_oneLine(lines);
    string strs_split = "self.listOfFinalCoefIndexes_arr = ";
    vector<vector<int>> listOfFinalCoefIndexes_arr = convert_str_dic_int(converted_line, strs_split);
    return listOfFinalCoefIndexes_arr;
}

vector<vector<vector<LagrangePolynomial>>> TuckerApproximation::getListOfBasicFcts( \
                                        vector<vector<double>> listOfDomainBorders, \
                                        vector<vector<double>>  listOfTuckerGridNodes, \
                                        vector<vector<vector<double>>> finalOrthNormalizedEigVects)
{
    vector<vector<vector<LagrangePolynomial>>> listOfBasicFctsUsingLagrangeInterpolation;
    vector<vector<LagrangePolynomial>> basisFcts_Axis_k;
    for (int Axis_k=0; Axis_k<dimension; Axis_k++)
    {
        basisFcts_Axis_k = getListOfInterpolationFcts(Axis_k, listOfDomainBorders, listOfTuckerGridNodes, \
          finalOrthNormalizedEigVects);
        listOfBasicFctsUsingLagrangeInterpolation[Axis_k] = basisFcts_Axis_k;
    }
    return listOfBasicFctsUsingLagrangeInterpolation;
}

/*
Les fichiers à lire :

GeneralInfor_UOXGd.txt (ou GeneralInfor_UOX.txt, GeneralInfor_MOX.txt)

TuckerInfor_macro_totale1.txt
TuckerInfor_macro_totale0.txt
TuckerInfor_macro_scattering011.txt
TuckerInfor_macro_scattering010.txt
TuckerInfor_macro_scattering001.txt
TuckerInfor_macro_scattering000.txt
TuckerInfor_macro_nu_fission1.txt
TuckerInfor_macro_nu_fission0.txt
TuckerInfor_macro_fission1.txt
TuckerInfor_macro_fission0.txt
TuckerInfor_macro_absorption1.txt
TuckerInfor_macro_absorption0.txt

*/
