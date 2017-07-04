#include "../include/Functions.hpp"

vector<string> Functions::allCoreTypes;
vector<string> Functions::allCrossSectionType;

Functions::Functions(string c, vector<string> cs)
{
    m_coreType = c;
    m_n = cs.size();
    m_crossSections = cs;
    m_tuckerProgram = make_shared<TuckerApproximation>(m_coreType,m_crossSections);
    m_directory = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_" + c + "_test";
    //cout << " - Directory of Tucker code : " << m_directory << endl;
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

void Functions::setCrossSectionType(vector<string> cs)
{
    m_crossSections = cs;
    m_n = m_crossSections.size();
}

void Functions::setAllCrossSectionType()
{
    for (string r : allCrossSectionType)
        m_crossSections.push_back(r);
    m_n = m_crossSections.size();
}

void Functions::setTuckerProgram()
{
    m_tuckerProgram = make_shared<TuckerApproximation>(m_coreType,m_crossSections);
}

vector<double> convert_str_to_vec(string s)
{
    s = Utils::replace(s, "(", "");
    s = Utils::replace(s, ")", "");
    s = Utils::replace(s, "[", "");
    s = Utils::replace(s, "]", "");
    s = Utils::replace(s, " ", "");
    vector<double> res;
    stringstream ss(s);
    string subs;
    size_t found;
    while (getline(ss, subs, ','))
    {
        found = subs.find("+");
        if (found!=string::npos)
            subs = subs.substr(0,found-1);
        res.push_back(stod(subs));
    }
    return res;
}

vector<double> parse_output(char* cmd)
{
    vector<double> res;
    char buf[BUFSIZE];
    FILE *fp;
    if ((fp = popen(cmd, "r")) == NULL)
        cerr << "Error opening pipe!" << endl;

    while (fgets(buf, BUFSIZE, fp) != NULL)
    {
        string buf_str(buf);
        cout << buf << endl;
        res = convert_str_to_vec(buf_str);
    }
    if (pclose(fp))
        cerr << "Command not found or exited with error status" << endl;
    return res;
}

vector<double> Functions::evaluate(MultiVariatePoint<double> x)
{

    //auto t0 = chrono::steady_clock::now();
    /*
    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";
    string cmd_s = "cd " + m_directory + "\npython eval.py ";
    for (int i=0; i<m_n; i++)
        cmd_s += m_crossSections[i] + " ";
    cmd_s += xstr;
    char *cmd = new char[cmd_s.length() + 1];
    strcpy(cmd, cmd_s.c_str());
    vector<double> res = parse_output(cmd);
    */
    //auto t1 = chrono::steady_clock::now();

    vector<double> result;
    for (string csName : m_crossSections)
        result.push_back(m_tuckerProgram->evaluate(x,csName));

    //auto t2 = chrono::steady_clock::now();

    //std::chrono::duration<double> delta1 = t1 - t0;
    //std::chrono::duration<double> delta2 = t2 - t1;
    /*
    cout << x << endl;
    cout << "[ ";
    for (size_t i=0; i<res.size()-1; ++i)
        cout << res[i] << " ; ";
    cout << res[res.size()-1] << " ]" << endl;
    cout << "[ ";
    for (size_t i=0; i<result.size()-1; ++i)
        cout << result[i] << " ; ";
    cout << result[result.size()-1] << " ]" << endl << endl;
    */
    return result;
}

bool Functions::validCoreType(string c)
{
    for (string s : allCoreTypes)
        if (s.compare(c)==0)
            return true;
    return false;
}

bool Functions::validCrossSections(string r)
{
  for (string s : allCrossSectionType)
      if (s.compare(r)==0)
          return true;
  return false;
}
