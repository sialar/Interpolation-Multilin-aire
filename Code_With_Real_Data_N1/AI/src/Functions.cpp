#include "../include/Functions.hpp"

vector<string> Functions::allCoreTypes;
vector<string> Functions::allCrossSectionType;

Functions::Functions(string c, string cs)
{
    m_coreType = c;
    m_crossSection = cs;
    m_directory = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_" + c + "_test";
}

void Functions::createFunctionsDataBase()
{
    allCoreTypes.push_back("MOX");
    allCoreTypes.push_back("UOX");
    allCoreTypes.push_back("UOX-Gd");

    allCrossSectionType.push_back("macro_absorption0");
    allCrossSectionType.push_back("macro_absorption1");
    allCrossSectionType.push_back("macro_fission0");
    allCrossSectionType.push_back("macro_fission1");
    allCrossSectionType.push_back("macro_nu*fission0");
    allCrossSectionType.push_back("macro_nu*fission1");
    allCrossSectionType.push_back("macro_scattering000");
    allCrossSectionType.push_back("macro_scattering001");
    allCrossSectionType.push_back("macro_scattering010");
    allCrossSectionType.push_back("macro_scattering011");
    allCrossSectionType.push_back("macro_totale0");
    allCrossSectionType.push_back("macro_totale1");
}

void Functions::setCoreType(string c)
{
    m_coreType = c;
    m_directory = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_" + c + "_test";
}

void Functions::setCrossSectionType(string cs)
{
    m_crossSection = cs;
}

double extractValue(string s)
{
    s = Utils::replace(s, "(", "");
    s = Utils::replace(s, ")", "");
    s = Utils::replace(s, " ", "");
    size_t found = s.find("+");
    if (found!=string::npos)
        s = s.substr(0,found-1);
    return stod(s);
}

double parse_output(char* cmd)
{
    double res = 0.0;
    char buf[BUFSIZE];
    FILE *fp;
    if ((fp = popen(cmd, "r")) == NULL)
        cerr << "Error opening pipe!" << endl;

    while (fgets(buf, BUFSIZE, fp) != NULL)
    {
        string buf_str(buf);
        res = extractValue(buf_str);
    }
    if (pclose(fp))
        cerr << "Command not found or exited with error status" << endl;
    return res;
}

double Functions::evaluate(MultiVariatePoint<double> x)
{
    string xstr = " ";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";
    string cmd_s = "cd " + m_directory + "\npython eval.py " + m_crossSection + xstr;
    char *cmd = new char[cmd_s.length() + 1];
    strcpy(cmd, cmd_s.c_str());
    return parse_output(cmd);
}

bool Functions::validCoreType(string c)
{
    for (string s : allCoreTypes)
        if (s.compare(c)==0)
            return true;
    return false;
}

bool Functions::validCrossSection(string r)
{
  for (string s : allCrossSectionType)
      if (s.compare(r)==0)
          return true;
  return false;
}
