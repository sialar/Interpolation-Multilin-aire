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


void compare(vector<double> u, vector<double> v, double threshold)
{
    int nb = 0;
    for (size_t i=0; i<u.size(); i++)
        if (abs(u[i]-v[i])>threshold)
        {
            nb++;
            cout << (i+1) << " " << u[i] << " " << v[i] << endl;
        }
    cout << nb << " diffrences found!" << endl;
}

vector<double> getFromFile(string fileName, int c)
{
  ifstream file_in(fileName, ios::out);
  vector<double> v, data;
  if(file_in)
  {
      string line;
      while (getline(file_in, line))
      {
          data = Utils::str2vector(line);
          v.push_back(data[c]);
      }
      file_in.close();
  }
  else cerr << "Error while opening the file!" << endl;
  return v;
}

int main( int argc, char* argv[] )
{
    string csName = "macro_scattering001";
    vector<double> u = getFromFile(Utils::projectPath + "fileName",0);
    vector<double> v = getFromFile(Utils::projectPath + "fileName_cpp",0);
    vector<double> w = getFromFile(Utils::projectPath + "fileName",0);

    compare(u, v, 1e-7);
    //compare(u, w, 1e-7);
    //compare(v, w, 1e-7);

    return 0;
}
