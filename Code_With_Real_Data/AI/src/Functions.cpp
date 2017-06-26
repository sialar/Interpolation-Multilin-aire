#include "../include/Functions.hpp"


Functions::Functions(int n, CoreType c, vector<ReactionType> vr)
{
    m_n = n;
    m_coreType = c;
    m_reactionTypes = vr;
    setFunctionsMap();
}

void Functions::setFunctionsMap()
{
    for (int c=MOX; c!=UOX_Gd; c++)
        for (int r=ABS_0; r!=TOT_1; r++)
            cout << c << "-" << r << endl;
}

void Functions::setReactionTypes(vector<ReactionType> vr)
{
    m_reactionTypes = vr;
    m_n = m_reactionTypes.size();
}


string Functions::getFunctionFileName(CoreType c, ReactionType r)
{
    string path = "/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/AI/data/";

    if (c==MOX) path += "mox/";
    else if (c==UOX) path += "uox/";
    else if (c==UOX_Gd) path += "uox-gd/";

    if (r==ABS_0) path += "abs0";
    else if (r==ABS_0) path += "abs0";
    else if (r==ABS_1) path += "abs1";
    else if (r==FIS_0) path += "fis0";
    else if (r==FIS_1) path += "fis1";
    else if (r==NU_FIS_0) path += "nu_fis0";
    else if (r==NU_FIS_1) path += "nu_fis1";
    else if (r==SCAT_00) path += "scat00";
    else if (r==SCAT_01) path += "scat01";
    else if (r==SCAT_10) path += "scat10";
    else if (r==SCAT_11) path += "scat11";
    else if (r==TOT_0) path += "tot0";
    else if (r==TOT_1) path += "tot1";
    return path;
}

vector<double> Functions::evaluate(MultiVariatePoint<double> x)
{
    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";
    system(("cd " + Utils::projectPath + "../Tucker/LIVRAISON_THESE_Paris6_UOX_test\n" \
            "python eval_UOX.py macro_totale0 " + xstr).c_str());
    string fileName = "/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/AI/data/cross_section_value.dat";
    ifstream file(fileName, ios::in);
    vector<double> res(1,0);
    if(file)
    {
        string line;
        getline(file,line);
        res[0] = stof(line);
        file.close();
    }
    else
        cerr << "Error while opening the file! " << endl;
    return res;
}

void Functions::updateValues()
{
  string fileName;

  for (int i=0; i<m_n; i++)
  {
      fileName = getFunctionFileName(m_coreType,m_reactionTypes[i]);
      ifstream file(fileName, ios::in);
      string line;
      if(file)
      {
          while (getline(file,line))
          {
          }
          file.close();
      }
      else
          cerr << "Error while opening the file! " << endl;

  }
}
