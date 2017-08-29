#include "../../include/Tucker/TuckerApproximation.hpp"

int TuckerApproximation::dimension = 5;
vector<string> TuckerApproximation::keys;

/******************************************************************************/
/******************************  utile functions ********************************/
/******************************************************************************/


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
string check_string(string strs_begin, string strs_end, string NameFile)
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
vector<string> get_list_check_string(string strs_begin, string strs_end, string NameFile)
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
string convert_multiLines_oneLine(string lines)
{
  return Utils::replace(lines,"\n"," ");
}
string str_key(string k)
{
  return "'" + k + "'";
}
vector<double> get_first_vect_from_str(string s)
{
    // s : " 'k' : [ a0, a1, a2, ..., an ], ... " ---> vec = map['k'] = [ a0, a1, a2, ..., an ]
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
vector<vector<double>> str_to_vect_of_vect(string s)
{
    vector<vector<double>> dic;
    s.pop_back();
    s.erase(0,1);

    string str_temp, strs = s;
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
         vector<double> vec;
         int i = 0;
         while (getline(ss, subs, ' '))
         {
              if (subs!="")
              {
                  //subs = Utils::getRealPart(subs);
                  vec.push_back(stod(subs));
                  i++;
              }
         }
         dic.push_back(vec);
    }
    return dic;
}
vector<double> str_to_vect(string s)
{
    s.pop_back();
    s.pop_back();
    s.erase(0,1);
    vector<double> dic;
    stringstream ss(s);
    string subs;
    while (getline(ss, subs, ' '))
    {
        if (subs!="") dic.push_back(stod(subs));
    }
    return dic;
}
vector<vector<int>> str_to_vect_of_vect_of_indexes(string s)
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
              vec.push_back(stoi(subs));
         dic.push_back(vec);
    }
    return dic;
}

map<string,vector<double>> str_to_vecs_map(string strs, string strs_split)
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
      dic.insert(pair<string,vector<double>>(k,get_first_vect_from_str(temp)));
    }
    return dic;
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<vector<double>> str_to_eigen_vects(string strs, string strs_split)
{
  vector<vector<double>> dic;
  vector<double> data;

  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = str_to_vect_of_vect(strs_dic);
    return dic;
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<double> str_to_tucker_coefs(string strs, string strs_split)
{
  vector<double> dic;
  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = str_to_vect(strs_dic);
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}
vector<vector<int>> str_to_coef_indexes(string strs, string strs_split)
{
  vector<vector<int>> dic;
  if (isSubString(strs,strs_split))
  {
    size_t found = strs.find(strs_split);
    string strs_dic = strs.substr(found+strs_split.length(),strs.length());
    dic = str_to_vect_of_vect_of_indexes(strs_dic);
  }
  else cerr << "Warning : strs_split is not in strs!" << endl;
  return dic;
}

/******************************************************************************/
/**************************  constuctor & methods  ****************************/
/******************************************************************************/
TuckerApproximation::TuckerApproximation(string _core, vector<string> listOfCrossSection)
{
    infoFileName = "Tucker/With_Real_Data/LIVRAISON_THESE_Paris6_" + _core + "_test/GeneralInfor_" + _core + ".txt";
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

TuckerApproximation::TuckerApproximation(string _core, string csName)
{
    infoFileName = "Tucker/With_Real_Data/LIVRAISON_THESE_Paris6_" + _core + "_test/GeneralInfor_" + _core + ".txt";
    core = _core;
    dimension = 5;
    keys.push_back("0");
    keys.push_back("1");
    keys.push_back("2");
    keys.push_back("3");
    keys.push_back("4");
    vector<string> listOfCrossSection(1,csName);
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
        listOfFiles.insert(pair<string,string>(cs,"Tucker/With_Real_Data/LIVRAISON_THESE_Paris6_" + core + "_test/TuckerInfor_" + s + ".txt"));
    }
}
void TuckerApproximation::setDomainBorders()
{
    string strs_begin = "listOfDomainBorders = {";
    string strs_end = "}";
    string lines = check_string(strs_begin, strs_end, infoFileName);
    string converted_line = convert_multiLines_oneLine(lines);
    converted_line = Utils::eraseExtraSpaces(converted_line);
    string strs_split = " = ";
    listOfDomainBorders = str_to_vecs_map(converted_line, strs_split);
}
void TuckerApproximation::setTuckerGridNodes()
{
    string strs_begin = "listOfTuckerGridNodes = {";
    string strs_end = "}";
    string lines = check_string(strs_begin, strs_end, infoFileName);
    string converted_line = convert_multiLines_oneLine(lines);
    converted_line = Utils::eraseExtraSpaces(converted_line);
    string strs_split = " = ";
    listOfTuckerGridNodes = str_to_vecs_map( converted_line, strs_split);
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
        string strs_split = "self.finalOrthNormalizedEigVects_Axis_k: ";
        orthNormalizedEigVects.insert(pair<string,vector<vector<double>>>(to_string(i), \
                                                  str_to_eigen_vects(converted_line, strs_split)));
    }
    return orthNormalizedEigVects;
}
vector<double> TuckerApproximation::getFinalTuckerCoeffs(string NameFile)
{
    string strs_begin = "Coefficient values";
    string strs_end = "]";
    string lines = check_string(strs_begin, strs_end, NameFile);

    string converted_line = convert_multiLines_oneLine(lines);
    converted_line = Utils::eraseExtraSpaces(converted_line);
    string strs_split = "self.FinalTuckerCoeffs =";
    vector<double> FinalTuckerCoeffs = str_to_tucker_coefs(converted_line, strs_split);
    return FinalTuckerCoeffs;
}
vector<vector<int>> TuckerApproximation::getListOfFinalCoefIndexes_arr(string NameFile)
{
    string strs_begin = "self.listOfFinalCoefIndexes_arr";
    string strs_end = "]]";
    string lines = check_string(strs_begin, strs_end, NameFile);
    string converted_line = convert_multiLines_oneLine(lines);
    converted_line = Utils::eraseExtraSpaces(converted_line);
    string strs_split = "self.listOfFinalCoefIndexes_arr = ";
    vector<vector<int>> listOfFinalCoefIndexes_arr = str_to_coef_indexes(converted_line, strs_split);
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
double TuckerApproximation::getInterpolation(vector<LagrangePolynomial> fct, double x)
{
  int n = fct.size();
  double min = Utils::min_elt(fct[0].axisValues());
  double max = Utils::max_elt(fct[n-1].axisValues());

  if (n == 1)
  {
    return fct[0].getInterpolation(x);
  }
  if (x < min) return fct[0].getInterpolation(x);
  if (x > max) return fct[n-1].getInterpolation(x);

  if (n > 1 && min <= x && x <= max)
  {
    int j = 0;
    while (j <= n-1)
    {
      double _min = Utils::min_elt(fct[j].axisValues());
      double _max = Utils::max_elt(fct[j].axisValues());
      if (_min <= x && x <= _max)
      {
        return fct[j].getInterpolation(x);
      }
      else j++;
    }
  }
  return x;
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
