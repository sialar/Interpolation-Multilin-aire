#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>
#include "../include/Utils.hpp"
#include "../include/Functions.hpp"
#include "../include/Tucker/LagrangePolynomial.hpp"
#include "../include/Tucker/TuckerApproximation.hpp"


using namespace std;


double convert_str_to_double(string s)
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
    return res[0];
}

double parse(char* cmd)
{
    double res = 0;
    char buf[BUFSIZE];
    FILE *fp;
    if ((fp = popen(cmd, "r")) == NULL)
        cerr << "Error opening pipe!" << endl;

    while (fgets(buf, BUFSIZE, fp) != NULL)
    {
        string buf_str(buf);
        res = convert_str_to_double(buf_str);
    }
    if (pclose(fp))
        cerr << "Command not found or exited with error status" << endl;
    return res;
}

double evaluate_py(MultiVariatePoint<double> x)
{
    string xstr = "";
    for (int i=0; i<x.getD(); i++)
        xstr += to_string(x(i)) + " ";
    string cmd_s = "cd /home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/Tucker/LIVRAISON_THESE_Paris6_MOX_test\npython eval.py macro_nu*fission1 ";
    cmd_s += xstr;
    char *cmd = new char[cmd_s.length() + 1];
    strcpy(cmd, cmd_s.c_str());
    return parse(cmd);
}

int main( int argc, char* argv[] )
{
    vector<string> csNames;
    csNames.push_back("macro_nu*fission1");
    TuckerApproximationPtr tucker = make_shared<TuckerApproximation>("MOX", csNames);

    MultiVariatePoint<double> t(5,0,0);
    t(0) = 0;
    t(1) = 20;
    t(2) = 0.40000000596;
    t(3) = 0;
    t(4) = 9.999999974752427e-07;

    vector<MultiVariatePoint<double>> points;
    vector<double> data, res, res_py;
    string csName = csNames[0];
    string s = Utils::replace(csName,"*","_");
    ifstream file("/home/sialar/Stage/LaboJ_LLions/Code/Code_With_Real_Data/Tucker/LIVRAISON_THESE_Paris6_MOX_test/FinalResults/" + s, ios::in);
    if(file)
    {
        string line;
        while (getline(file, line))
        {
            MultiVariatePoint<double> p(5,0,0);
            data = Utils::str2vector(line);
            for (int i=0; i<5; i++)
                p(i) = data[i+1];
            res.push_back(data[8]);
            points.push_back(p);
        }
        file.close();
    }
    else cerr << "Error while opening " << csName << " file!" << endl;

    for (size_t i=2000; i<2100; i++)
    {
        double x = tucker->evaluate(points[i],"macro_nu*fission1");
        double y = res[i];
        double z = evaluate_py(points[i]);
        cout << i << " " << points[i] << " " << x << " " << y << " " << z <<endl;
    }
    return 0;
}
