#include "../include/Interpolation.hpp"

template <typename T>
Interpolation<T>::Interpolation(int d, int nIter) : m_d(d), m_maxIteration(nIter)
{
    m_interpolationPoints.resize(m_d);
}

template <typename T>
void Interpolation<T>::displayInterpolationPointsInEachDirection()
{
    vector<double>::iterator it;
    for (int i=0; i<m_d; ++i)
    {
        cout << " - " << m_interpolationPoints[i].size() << " points in direction " << i << " : { ";
        for (it=m_interpolationPoints[i].begin(); it!=m_interpolationPoints[i].end(); it++)
            cout << *it << " ";
        cout << "}" << endl;
    }
}

template <typename T>
void Interpolation<T>::displayInterpolationMultiVariatePoints()
{
    cout << " - Interpolation nodes: { ";
    for (MultiVariatePoint<double> x : m_interpolationNodes)
        cout << x << " ";
    cout << "}" << endl;
}

template <typename T>
void Interpolation<T>::saveErrorsInFile()
{
  ofstream file("data/interpolation_error.txt", ios::out | ios::trunc);
  if(file)
  {
      file << m_errors.size() << endl;
      map<int,double>::iterator it;
      for (it=m_errors.begin(); it!=m_errors.end(); it++)
          file << get<0>(*it) << " " << get<1>(*it) << endl;
      file.close();
  }
  else
      cerr << "Error while opening the file!" << endl;
}
