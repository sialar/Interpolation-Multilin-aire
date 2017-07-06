#include "../include/Functions.hpp"

vector<vector<vector<double>>> Functions::m_coefs;
int Functions::m_polynomialDegree;
int Functions::m_nbPerturbations;
vector<double> Functions::m_perturbations;
vector<double> Functions::m_perturbationsWidth;

string replace(string strs, string str_old, string str_new)
{
  size_t found = strs.find(str_old);
  while (found!=string::npos)
  {
    strs.replace(found, str_old.length(), str_new);
    found = strs.find(str_old);
  }
  return strs;
}
vector<double> convert_str_to_vec(string s)
{
    s = replace(s, "(", "");
    s = replace(s, ")", "");
    s = replace(s, "[", "");
    s = replace(s, "]", "");
    s = replace(s, " ", "");
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
        res = convert_str_to_vec(buf_str);
    }
    if (pclose(fp))
        cerr << "Command not found or exited with error status" << endl;
    return res;
}

double convertToFunctionDomain(double a, double b, double x)
{
  return x*(b-a)/2 + (a+b)/2;
}

vector<double> Functions::evaluate(MultiVariatePoint<double> x, int n)
{
    vector<vector<double>> dom(5);
    vector<double> interval(2);
    interval[0] = 0;
    interval[1] = 80000;
    dom[0] = interval;
    interval[0] = 20;
    interval[1] = 1800;
    dom[1] = interval;
    interval[0] = 0.4;
    interval[1] = 1;
    dom[2] = interval;
    interval[0] = 0;
    interval[1] = 1.74835e-05;
    dom[3] = interval;
    interval[0] = 1e-06;
    interval[1] = 2;
    dom[4] = interval;

    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(convertToFunctionDomain(dom[i][0],dom[i][1],x(i))) + " ";
    string cmd_s = "cd /home/sialar/Stage/LaboJ_LLions/Code/Code_With_Analytic_Function/MOX\npython eval.py macro_nu*fission0 ";
    cmd_s += xstr;
    char *cmd = new char[cmd_s.length() + 1];
    strcpy(cmd, cmd_s.c_str());
    vector<double> res = parse_output(cmd);
    //for (int i=0; i<x.getD(); i++)
    //    cout << x(i) << " " << to_string(convertToFunctionDomain(dom[i][0],dom[i][1],x(i))) << " " << convertToFunctionDomain(0,80000,-1) << endl;


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
    return res;
}

double Functions::getPointInPerturbationNeighborhood()
{
    if (!m_perturbations.size()) return 0;
    return m_perturbations[0] + m_perturbationsWidth[0]/10;
}

double Functions::hat(double a, double b, double fa, double fb, double t)
{
    double c = (a+b)/2;
    if (t <= c && t >= a) return ((t-a)/(c-a)+fa);
    else if (t >= c && t <= b) return ((t-b)/(c-b)+fb);
    else return 0;
}

double Functions::changeFunctionDomain(double a, double b, double x)
{
    return (2*x)/(b-a) + 1 - (2*b)/(b-a);
}

vector<double> Functions::f(MultiVariatePoint<double> x, int n)
{
    // write here the interpolated function
    return autoPolynomialFunction(x,n);
}

vector<double> Functions::toAlternatingVector(double x, int n)
{
    vector<double> res(n,x);
    for (int i=0; i<n; i++)
        if (i%2) res[i] *= -1;
    return res;
}

vector<double> Functions::functionToPlot(MultiVariatePoint<double> x, int n)
{
    double temp = 1;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    return toAlternatingVector(sqrt(1-x(0)*x(0)) * temp, n);
}

vector<double> Functions::h(MultiVariatePoint<double> x, int n)
{
    double temp = 0;
    for (int i=0; i<x.getD(); i++) temp += x(i);
    if (x.getD()>1) return toAlternatingVector(sin(temp),n);
    else return toAlternatingVector(2*exp(-x(0)*x(0) + sin(2*M_PI*x(0))),n);
}

int Functions::inOneHat(double x)
{
    double xc, xe;
    for (int i=0; i<m_nbPerturbations; i++)
    {
        xc = m_perturbations[i];
        xe = m_perturbationsWidth[i];
        if (x>=xc-xe/2 && x<=xc+xe/2)
            return i;
    }
    return -1;
}
vector<double> Functions::phi(MultiVariatePoint<double> x, int n)
{
    double temp = 1, f = 1;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    int hatPos = inOneHat(x(0));
    double xc, xe, a, b;
    if (hatPos>-1)
    {
        xc = m_perturbations[hatPos];
        xe = m_perturbationsWidth[hatPos];
        a = xc - xe/2;
        b = xc + xe/2;
        return toAlternatingVector(hat(a,b,sin(2*M_PI*f*a),sin(2*M_PI*f*b),x(0)) * temp,n);
    }
    else
        return toAlternatingVector(sin(2*M_PI*f*x(0)) * temp,n);
}

vector<double> Functions::cosinus(MultiVariatePoint<double> x, int n)
{
    double temp = 1, f = 2;
    for (int i=1; i<x.getD(); i++) temp *= exp(-x(i));
    return toAlternatingVector(cos(2*M_PI*f*x(0)) * temp,n);
}

vector<double> Functions::function1D(MultiVariatePoint<double> x, int n)
{
    if (x.getD()==1)
    {
        double f = 10;
        if (x(0)<0) return toAlternatingVector(cos(2*M_PI*f*x(0)),n);
        else return toAlternatingVector(cos(0.1*M_PI*f*x(0)),n);
    }
    else return toAlternatingVector(0,n);
}

vector<double> Functions::sinOfNorm2(MultiVariatePoint<double> x, int n)
{
    double temp = 0;
    for (int i=0; i<x.getD(); i++)
        temp += pow(x(i),2);
    return toAlternatingVector(sin(sqrt(temp)),n);
}

vector<double> Functions::autoPolynomialFunction(MultiVariatePoint<double> x, int n)
{
    vector<double> res;
    res.resize(n,0);
    for (int i=0; i<n; i++)
        for (int j=0; j<m_polynomialDegree; j++)
            for (int k=0; k<x.getD(); k++)
                res[i] += m_coefs[i][j][k] * pow(x(k),j+1);
    return res;
}

void Functions::setCoefs(int degree, int d, int n)
{
    m_nbPerturbations = Utils::randomValue(1,3);
    for (int i=0; i<m_nbPerturbations; i++)
    {
        m_perturbations.push_back(Utils::randomValue(0.1,0.9));
        m_perturbationsWidth.push_back(Utils::randomValue(0.01,0.05));
    }
    m_polynomialDegree = degree;
    m_coefs.resize(n);
    for (int i=0; i<n; i++)
        m_coefs[i].resize(degree);
    for (int i=0; i<n; i++)
        for (int j=0; j<degree; j++)
            m_coefs[i][j].resize(d);

    for (int i=0; i<n; i++)
        for (int j=0; j<degree; j++)
            for (int k=0; k<d; k++)
                m_coefs[i][j][k] = Utils::randomValue(-1,1);

    saveCoefsInFile(d,n);
}

void Functions::saveCoefsInFile(int d, int n)
{
    ofstream file(Utils::projectPath + "data/polynomial_coefs.txt", ios::out | ios::trunc);
    if(file)
    {
        file << n << " " << m_polynomialDegree << " " << d << endl;
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<m_polynomialDegree; ++j)
            {
                for (int k=0; k<d; k++)
                {
                    file << m_coefs[i][j][k] << " ";
                }
                file << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}

void Functions::savePerturbationsInFile()
{
    ofstream file(Utils::projectPath + "data/perturbations.txt", ios::out | ios::trunc);
    if(file)
    {
        file << m_nbPerturbations << endl;

        for (int i=0; i<m_nbPerturbations; ++i)
            file << m_perturbations[i] << " ";
        file << endl;

        for (int i=0; i<m_nbPerturbations; ++i)
            file << m_perturbationsWidth[i] << " ";
        file << endl;

        file.close();
    }
    else
        cerr << "Erreur à l'ouverture du fichier!" << endl;
}
