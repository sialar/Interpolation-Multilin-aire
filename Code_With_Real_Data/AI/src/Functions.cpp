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

void Functions::setTuckerProgram()
{
    m_tuckerProgram = make_shared<TuckerApproximation>(m_coreType,m_reactionTypes);
}

vector<double> convert_str_to_vec(string s)
{
    s.pop_back();
    s.erase(0,1);
    vector<double> res;
    stringstream ss(s);
    string subs;
    size_t found;
    while (getline(ss, subs, ','))
    {
        found = subs.find("+");
        if (found!=string::npos)
        {
            subs = subs.substr(0,found-1);
            subs.erase(0,1);
        }
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
        res = convert_str_to_vec(buf_str);
    }
    if (pclose(fp))
        cerr << "Command not found or exited with error status" << endl;
    return res;
}

vector<double> Functions::evaluate(MultiVariatePoint<double> x)
{
    /*
    auto t0 = chrono::steady_clock::now();

    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";
    string cmd_s = "cd " + m_directory + "\npython eval.py ";
    for (int i=0; i<m_n; i++)
        cmd_s += m_reactionTypes[i] + " ";
    cmd_s += xstr;
    char *cmd = new char[cmd_s.length() + 1];
    strcpy(cmd, cmd_s.c_str());
    vector<double> res = parse_output(cmd);

    auto t1 = chrono::steady_clock::now();
    */
    vector<double> result;
    for (string csName : m_reactionTypes)
        result.push_back(m_tuckerProgram->evaluate(x,csName));
    /*
    auto t2 = chrono::steady_clock::now();

    std::chrono::duration<double> delta1 = t1 - t0;
    std::chrono::duration<double> delta2 = t2 - t1;
    cout << x << endl;
    cout << "[ ";
    for (size_t i=0; i<res.size()-1; ++i)
        cout << m_reactionTypes[i] << "|" << res[i] << " ; ";
    cout << res[res.size()-1] << " ]" << delta1.count() << endl;
    cout << "[ ";
    for (size_t i=0; i<result.size()-1; ++i)
        cout << m_reactionTypes[i] << "|" << result[i] << " ; ";
    cout << result[result.size()-1] << " ]" << delta2.count() << endl << endl;
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

bool Functions::validReactionType(string r)
{
  for (string s : allReactionTypes)
      if (s.compare(r)==0)
          return true;
  return false;
}
