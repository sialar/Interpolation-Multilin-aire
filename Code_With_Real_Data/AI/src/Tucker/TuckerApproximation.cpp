#include "../../include/Tucker/TuckerApproximation.hpp"

int TuckerApproximation::dimension = 5;
vector<string> TuckerApproximation::keys;

TuckerApproximation::TuckerApproximation(string _core, vector<string> listOfCrossSection)
{
    infoFileName = Utils::projectPath + "AI/data/" + _core + "/GeneralInfor_" + _core + ".txt";
    core = _core;
    dimension = 5;
    keys.push_back("0");
    keys.push_back("1");
    keys.push_back("2");
    keys.push_back("3");
    keys.push_back("4");
    setListOfCrossSection(listOfCrossSection);
    setListOfFiles();
    setDomainBorders();
    setTuckerGridNodes();
    setAll();
}
void TuckerApproximation::setListOfCrossSection(vector<string> listOfCrossSection)
{
    listOfCrossSectionNames = listOfCrossSection;
}
void TuckerApproximation::setListOfFiles()
{
    for (string cs : listOfCrossSectionNames)
    {
        string s = Utils::replace(cs,"*","_");
        listOfFiles.insert(pair<string,string>(cs,Utils::projectPath + "AI/data/" + core + "/TuckerInfor_" + s + ".txt"));
    }
}
void TuckerApproximation::setDomainBorders()
{
    string strs_begin = "listOfDomainBorders = {";
    string strs_end = "}";
    string lines = check_string(strs_begin, strs_end, infoFileName);
    string converted_line = convert_multiLines_oneLine(lines);
    string strs_split = " = ";
    listOfDomainBorders = convert_str_dic_map(converted_line, strs_split);
}
void TuckerApproximation::setTuckerGridNodes()
{
    string strs_begin = "listOfTuckerGridNodes = {";
    string strs_end = "}";
    string lines = check_string(strs_begin, strs_end, infoFileName);
    string converted_line = convert_multiLines_oneLine(lines);
    string strs_split = " = ";
    listOfTuckerGridNodes = convert_str_dic_map( converted_line, strs_split);
}
void TuckerApproximation::setAll()
{
    for (string csName : listOfCrossSectionNames)
    {
        string NameFile = listOfFiles[csName];
        finalOrthNormalizedEigVects.insert(pair<string,map<string,vector<vector<double>>>>(csName,getFinalOrthNormalizedEigVects(NameFile)));
        FinalTuckerCoeffs.insert(pair<string,vector<double>>(csName, getFinalTuckerCoeffs(NameFile)));
        listOfFinalCoefIndexes_arr.insert(pair<string,vector<vector<int>>>(csName, getListOfFinalCoefIndexes_arr(NameFile)));
        listOfBasicFctsUsingLagrangeInterpolation.insert(pair<string,map<string,vector<vector<LagrangePolynomial>>>>(csName, getListOfBasicFcts(csName)));
    }
}

/******************************************************************************/
/*************************  Intermediate functions ****************************/
/******************************************************************************/

double min_elt(vector<double> v)
{
    double m = v[0];
    for (double x : v)
        if (x<m)
            m = x;
    return m;
}
double max_elt(vector<double> v)
{
    double m = v[0];
    for (double x : v)
        if (x>m)
            m = x;
    return m;
}
string str_key(string k)
{
  return "'" + k + "'";
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
bool isSubString(string str, string substr)
{
    size_t found = str.find(substr);
    return (found!=string::npos);
}
vector<double> getFirstVectorFromString(string s)
{
    vector<double> vec;
    size_t found = s.find("[");
    s.replace(0, found+1, "");
    found = s.find("]");
    s.replace(found, s.length(), "");
    stringstream ss(s);
    string subs;
    int i = 0;
    while (getline(ss, subs, ','))
    {
        vec.push_back(stod(subs));
        i++;
    }
    return vec;
}
vector<vector<double>> convertStringToVectorOfVectors(string s)
{
    s.pop_back();
    s.erase(0,1);
    string str_temp, strs = s;
    vector<vector<double>> dic;
    size_t foundj, found2, found1 = strs.find("[");
    while (found1!=string::npos)
    {
         strs = strs.substr(found1+1,strs.length());
         found2 = strs.find("]");
         str_temp = strs;
         str_temp = str_temp.substr(0,found2);
         found1 = strs.find("[");

         stringstream ss(str_temp);
         string subs;
         vector<double> vec;
         int i = 0;
         //getline(ss, subs, ',');
         while (getline(ss, subs, ','))
         {
              if (subs!="")
              {
                  foundj = strs.find("+");
                  if (foundj!=string::npos)
                      subs = subs.substr(0,foundj-1);
                  if (subs!="")  vec.push_back(stod(subs));
                  i++;
              }
         }
         dic.push_back(vec);
    }
    return dic;
}
vector<double> convertStringToVector(string s)
{
    s.pop_back();
    s.pop_back();
    s.erase(0,2);
    s = Utils::eraseExtraSpaces(s);
    s = Utils::replace(s, " ", ",");
    vector<double> dic;
    stringstream ss(s);
    string subs;
    while (getline(ss, subs, ','))
    {
        if (subs!="") dic.push_back(stod(subs));
    }
    return dic;
}
vector<vector<int>> convertStringToVectorOfVectorsOfIndexes(string s)
{
    s.pop_back();
    s.erase(0,1);
    string str_temp, strs = s;
    vector<vector<int>> dic;
    size_t found2, found1 = strs.find("[");
    while (found1!=string::npos)
    {
         strs = strs.substr(found1+1,strs.length());
         found2 = strs.find("]");
         str_temp = strs;
         str_temp = str_temp.substr(0,found2);
         found1 = strs.find("[");

         stringstream ss(str_temp);
         string subs;
         vector<int> vec;
         while (getline(ss, subs, ','))
              vec.push_back(stod(subs));
         dic.push_back(vec);
    }
    return dic;
}

/******************************************************************************/
/************************* static functions ***********************************/
/******************************************************************************/

string TuckerApproximation::check_string(string strs_begin, string strs_end, string NameFile)
{
  ifstream file(NameFile, ios::in);
  if(file)
  {
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
        file.close();
        return found_strs;
      }
    }
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
        while (!isSubString(line,strs_end))
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
map<string,vector<double>> TuckerApproximation::convert_str_dic_map(string strs, string strs_split)
{
  map<string,vector<double>> dic;
  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string temp, strs_dic = strs.substr(found+strs_split.length(),strs.length());
    for (string k : TuckerApproximation::keys)
    {
      found = strs_dic.find(str_key(k));
      temp = strs_dic.substr(found,strs_dic.length());
      dic.insert(pair<string,vector<double>>(k,getFirstVectorFromString(temp)));
    }
    return dic;
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<vector<double>> TuckerApproximation::convert_str_dic_eigVects(string strs, string strs_split)
{
  vector<vector<double>> dic;
  vector<double> data;

  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = convertStringToVectorOfVectors(strs_dic);
    return dic;
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<double> TuckerApproximation::convert_str_dic_tucker_coef(string strs, string strs_split)
{
  vector<double> dic;
  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = convertStringToVector(strs_dic);
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<vector<int>> TuckerApproximation::convert_str_dic_coef_index(string strs, string strs_split)
{
  vector<vector<int>> dic;
  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = convertStringToVectorOfVectorsOfIndexes(strs_dic);
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
double TuckerApproximation::getInterpolation(vector<LagrangePolynomial> fct, double x)
{
  int n = fct.size();
  double min = min_elt(fct[0].axisValues());
  double max = max_elt(fct[n-1].axisValues());

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
vector<double> TuckerApproximation::getInterpolationArr(vector<LagrangePolynomial> fct, MultiVariatePoint<double> x)
{
  vector<double> listOfValues;
  for (int i=0; i<x.getD(); i++) listOfValues.push_back(getInterpolation(fct,x(i)));
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
  double min_dist = numeric_limits<double>::max();
  double argmin = 0;
  for (double x : myList)
  {
    if (abs(x-value)<min_dist)
    {
      min_dist = abs(x-value);
      argmin = x;
    }
  }
  return argmin;
}

/******************************************************************************/
/************************* main functions *************************************/
/******************************************************************************/
double TuckerApproximation::evaluate(MultiVariatePoint<double> point, string csName)
{
  double approxValue = 0.0;
  for (int i=0; i<int(listOfFinalCoefIndexes_arr[csName].size()); i++)
  {
    double d = 1.0;
    vector<LagrangePolynomial> fct;
    for (int Axis_k=0; Axis_k<dimension; Axis_k++)
    {
      fct = listOfBasicFctsUsingLagrangeInterpolation[csName][to_string(Axis_k)][listOfFinalCoefIndexes_arr[csName][i][Axis_k]];
      d = d*getInterpolation(fct, point(Axis_k));
    }
    approxValue = approxValue +  FinalTuckerCoeffs[csName][i]*d;
  }
  return approxValue;
}
vector<vector<LagrangePolynomial>> TuckerApproximation::getListOfInterpolationFcts(string Axis_k, string csName)
{
    vector<vector<LagrangePolynomial>> listOfBasicFctsUsingLagrangeInterpolation_Axis_k;
    int nbExtremes =  listOfDomainBorders[Axis_k].size();
    int nbOfFcts_Axis_k = finalOrthNormalizedEigVects[csName][Axis_k].size();
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
                                          vectorSequence(finalOrthNormalizedEigVects[csName][Axis_k][i],j1,j2));

                polLagrange_ki.push_back(interp_k_couplage);
            }
        }
        else if  (nbExtremes == 2)
        {
            LagrangePolynomial interp(listOfTuckerGridNodes[Axis_k],finalOrthNormalizedEigVects[csName][Axis_k][i]);
            polLagrange_ki.push_back(interp);
        }
        listOfBasicFctsUsingLagrangeInterpolation_Axis_k.push_back(polLagrange_ki);
    }
    return listOfBasicFctsUsingLagrangeInterpolation_Axis_k;
}
map<string,vector<vector<double>>> TuckerApproximation::getFinalOrthNormalizedEigVects(string NameFile)
{
    string strs_begin = "Orthonormalized eigenvectors";
    string strs_end = "]]";
    vector<string> list_check_string = get_list_check_string(strs_begin, strs_end, NameFile);
    map<string,vector<vector<double>>> orthNormalizedEigVects;
    for (int i=0; i<dimension; i++)
    {
        string lines = list_check_string[i];
        string converted_line = convert_multiLines_oneLine(lines);
        converted_line = Utils::eraseExtraSpaces(converted_line);
        converted_line = Utils::replace(converted_line, " ", ",");
        string strs_split = "self.finalOrthNormalizedEigVects_Axis_k:,";
        orthNormalizedEigVects.insert(pair<string,vector<vector<double>>>(to_string(i), \
                                                  convert_str_dic_eigVects(converted_line, strs_split)));
    }
    return orthNormalizedEigVects;
}
vector<double> TuckerApproximation::getFinalTuckerCoeffs(string NameFile)
{
    string strs_begin = "Coefficient values";
    string strs_end = "]";
    string lines = check_string(strs_begin, strs_end, NameFile);

    string converted_line = convert_multiLines_oneLine(lines);
    string strs_split = "self.FinalTuckerCoeffs =";
    vector<double> FinalTuckerCoeffs = convert_str_dic_tucker_coef(converted_line, strs_split);
    return FinalTuckerCoeffs;
}
vector<vector<int>> TuckerApproximation::getListOfFinalCoefIndexes_arr(string NameFile)
{
    string strs_begin = "self.listOfFinalCoefIndexes_arr";
    string strs_end = "]]";
    string lines = check_string(strs_begin, strs_end, NameFile);
    string converted_line = convert_multiLines_oneLine(lines);
    string strs_split = "self.listOfFinalCoefIndexes_arr = ";
    vector<vector<int>> listOfFinalCoefIndexes_arr = convert_str_dic_coef_index(converted_line, strs_split);
    return listOfFinalCoefIndexes_arr;
}
map<string,vector<vector<LagrangePolynomial>>> TuckerApproximation::getListOfBasicFcts(string csName)
{
    map<string,vector<vector<LagrangePolynomial>>> listOfBasicFctsUsingLagrangeInterpolation;
    vector<vector<LagrangePolynomial>> basisFcts_Axis_k;
    for (string Axis_k : TuckerApproximation::keys)
    {
        basisFcts_Axis_k = getListOfInterpolationFcts(Axis_k,csName);
        listOfBasicFctsUsingLagrangeInterpolation.insert(pair<string,vector<vector<LagrangePolynomial>>>(Axis_k,basisFcts_Axis_k));
    }
    return listOfBasicFctsUsingLagrangeInterpolation;
}
