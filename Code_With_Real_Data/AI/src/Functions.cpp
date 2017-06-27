#include "../include/Functions.hpp"

vector<string> Functions::allCoreTypes;
vector<string> Functions::allReactionTypes;

Functions::Functions(int n, string c, vector<string> vr)
{
    m_n = n;
    m_coreType = c;
    m_reactionTypes = vr;
    m_directory = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_" + c + "_test";
    cout << " - Directory of Tucker code : " << m_directory << endl;
}

void Functions::createFunctionsDataBase()
{
    allCoreTypes.push_back("MOX");
    allCoreTypes.push_back("UOX");
    allCoreTypes.push_back("UOX-Gd");

    allReactionTypes.push_back("macro_absorption0");
    allReactionTypes.push_back("macro_absorption1");
    allReactionTypes.push_back("macro_fission0");
    allReactionTypes.push_back("macro_fission1");
    allReactionTypes.push_back("macro_nu*fission0");
    allReactionTypes.push_back("macro_nu*fission1");
    allReactionTypes.push_back("macro_scattering000");
    allReactionTypes.push_back("macro_scattering001");
    allReactionTypes.push_back("macro_scattering010");
    allReactionTypes.push_back("macro_scattering011");
    allReactionTypes.push_back("macro_totale0");
    allReactionTypes.push_back("macro_totale1");
}

void Functions::setCoreType(string c)
{
    m_coreType = c;
    m_directory = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_" + c + "_test";
    cout << " - Directory of Tucker code : " << m_directory << endl;
}

void Functions::setReactionTypes(vector<string> vr)
{
    m_reactionTypes = vr;
    m_n = m_reactionTypes.size();
}

void Functions::setAllReactionTypes()
{
    for (string r : allReactionTypes)
        m_reactionTypes.push_back(r);
    m_n = m_reactionTypes.size();
}

vector<double> Functions::evaluate(MultiVariatePoint<double> x)
{
    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";

    string fileName = Utils::projectPath + "AI/data/cross_section_value.dat";
    string cmd = "cd " + m_directory + "\npython eval.py ";
    for (int i=0; i<m_n; i++)
        cmd += m_reactionTypes[i] + " ";
    cmd += xstr;

    if (system(cmd.c_str()) == -1)
    {
        cout << " - Failed to execute Tucker program!" << endl;
        exit(0);
    }

    ifstream file(fileName, ios::in);
    vector<double> res;
    if(file)
    {
        string line;
        while (getline(file,line))
            res.push_back(stof(line));
        file.close();
    }
    else
        cerr << "Error while opening the file! " << endl;
    return res;
}

bool Functions::validCoreType(string c)
{
    for (string s : allCoreTypes)
        if (s.compare(c)==0)
            return true;
    return false;
}

bool Functions::validReactionType(string r)
{
  for (string s : allReactionTypes)
      if (s.compare(r)==0)
          return true;
  return false;
}
