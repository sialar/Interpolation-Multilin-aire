#include "../include/Dichotomy.hpp"

int Dichotomy::getIndice(double l)
{
    if (l == 0) return 0;
    if (l == -1)  return 1;
    if (l == 1)  return 2;
    double temp = l;
    int n = 0;
    while (floor(temp) != temp)
    {
        temp *= 2;
        n++;
    }
    int start = 2;
    for (int k=1; k<n; k++)
        start += pow(2,n-k);
    double a = 0.5 * (temp + pow(2,n)) + 1;
    return start + int(a);
}
double Dichotomy::getValue(int i)
{
  if (i == 0) return 0;
  if (i == 1)  return -1;
  if (i == 2)  return 1;

  int temp = i - 2; // nb points a partir de -0.5
  int e = 1;
  while (temp > pow(2,e))
  {
    temp -= pow(2,e);
    e++;
  }
  return (2*temp - 1 - pow(2,e)) / pow(2,e);
}
int Dichotomy::getParentIndice(int i)
{
    if (i < 3) return 0;
    if (i == 3) return 1;
    if (i == 4) return 2;

    int temp = i - 2;
    int e = 1;
    while (temp > pow(2,e))
    {
        temp -= pow(2,e);
        e++;
    }
    // le noeud i est le tempème noeud dans le niveau e
    // son parent est le numéro floor((temp+1)/2) dans le niveau (e-1)
    return (pow(2,e-1) + ((temp%2) ? (temp+1) : temp)/2);
}
double Dichotomy::getParentValue(int i)
{
    if (!i) return 2;
    return getValue(getParentIndice(i));
}
vector<double> Dichotomy::computeChildrenValue(int i)
{
    double parentVal, curentVal, offset;
    vector<double> childrenVal;
    parentVal = getParentValue(i);
    curentVal = getValue(i);
    offset = abs(parentVal - curentVal) / 2;
    if (i==1) childrenVal.push_back(curentVal + offset);
    else if (i==2) childrenVal.push_back(curentVal - offset);
    else
    {
        childrenVal.push_back(curentVal + offset);
        childrenVal.push_back(curentVal - offset);
    }
    return childrenVal;
}
bool Dichotomy::isAncestorOf(int i, int j) // i ancestor of j
{
    if (i>j) return false;
    if (i==j || i==0) return true;
    int temp = getParentIndice(j);
    while (temp!=0)
    {
        if (i==temp) return true;
        temp = getParentIndice(temp);
    }
    return false;
}
