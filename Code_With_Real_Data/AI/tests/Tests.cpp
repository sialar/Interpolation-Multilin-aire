#include <python3.5/Python.h>
#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;

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

int main( int argc, char* argv[] )
{
    string s = "Orthonormalized eigenvectors, self.finalOrthNormalizedEigVects_Axis_k: [[0.00348374+0.j  0.00348397+0.j  0.00348418+0.j  0.00348418+0.j   0.00348511+0.j  0.00348728+0.j  0.00348990+0.j  0.00349268+0.j   0.00349535+0.j  0.00349757+0.j  0.00349904+0.j  0.00349955+0.j   0.00349955+0.j  0.00350317+0.j  0.00351319+0.j  0.00352745+0.j   0.00354309+0.j  0.00355741+0.j  0.00356847+0.j  0.00357531+0.j   0.00357762+0.j] [0.00572681+0.j  0.00572708+0.j  0.00573095+0.j  0.00573095+0.j   0.00576145+0.j  0.00575093+0.j  0.00559647+0.j  0.00535024+0.j   0.00508314+0.j  0.00485008+0.j  0.00469311+0.j  0.00463768+0.j   0.00463768+0.j  0.00424026+0.j  0.00309074+0.j  0.00134258+0.j  -0.00073086+0.j -0.00280226+0.j -0.00454485+0.j -0.00569967+0.j  -0.00610270+0.j]]";
    string str_old = " ";
    string str_new = ",";
    string res = replace(s,str_old,str_new);
    cout << s << endl;
    cout << res << endl;
    return 0;
}
