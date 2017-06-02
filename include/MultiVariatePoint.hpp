#ifndef MULTIVARIATEPOINT
#define MULTIVARIATEPOINT

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <memory>
#include <limits>

#include "BinaryTree.hpp"

using namespace std;

template <typename T>
class MultiVariatePoint
{
    private:
        int m_d;
        T* m_nu = NULL;
        double m_alpha;
        double m_initialAlpha;
        int m_waitingTime;

   public:

        ~MultiVariatePoint();
        MultiVariatePoint() : m_d(0) {};
        MultiVariatePoint(int d, T val);
        MultiVariatePoint(const MultiVariatePoint<T>& nu);

        int getD() const { return m_d; };

        int getWaitingTime() { return m_waitingTime; };
        void incrWaitingTime() { m_waitingTime++; };

        double getAlpha() const { return m_alpha; };
        void setAlpha(double alpha) { m_alpha = alpha; };

        bool alphaAlreadyComputed() { return (m_alpha != m_initialAlpha); };
        void reinit() { m_alpha = m_initialAlpha; };
        
        static MultiVariatePoint<T> toMultiVariatePoint(vector<T> vec);
        static MultiVariatePoint<T> toBiVariatePoint(T vec0, T vec1);
        static MultiVariatePoint<T> toMonoVariatePoint(T vec);

        T &operator()(int d) const { return m_nu[d]; };
        MultiVariatePoint& operator+=(const MultiVariatePoint & v);
        MultiVariatePoint& operator=(MultiVariatePoint const &nu);
};

template <typename T>
using MultiVariatePointPtr = shared_ptr<MultiVariatePoint<T>>;

template <typename T>
MultiVariatePoint<T>::~MultiVariatePoint()
{
    if (m_nu != NULL)
        delete[] m_nu;
}

template <typename T>
MultiVariatePoint<T>::MultiVariatePoint(int d, T val) : m_d(d)
{
    if (m_d != 0)
    {
        m_nu = new T[m_d];
        for (int i=0; i<m_d; i++)
            m_nu[i] = val;
    }
    m_alpha = numeric_limits<int>::max();
    m_initialAlpha = numeric_limits<int>::max();
    m_waitingTime = 0;
}

template <typename T>
MultiVariatePoint<T>::MultiVariatePoint(const MultiVariatePoint<T>& nu)
{
    m_d = nu.m_d;
    if (m_d != 0)
    {
        m_nu = new T[m_d];
        for (int i=0; i<m_d; i++)
            m_nu[i] = nu.m_nu[i];
    }
    m_alpha = numeric_limits<int>::max();
    m_initialAlpha = numeric_limits<int>::max();
    m_waitingTime = 0;
};
template <typename T>
MultiVariatePoint<T> MultiVariatePoint<T>::toMultiVariatePoint(vector<T> vec)
{
    MultiVariatePoint<T> x(vec.size(),0);
    for (int i=0; i<x.getD(); i++)
       x(i) = vec[i];
    return x;
}
template <typename T>
MultiVariatePoint<T> MultiVariatePoint<T>::toBiVariatePoint(T t0, T t1)
{
    MultiVariatePoint<T> x(2,0);
    x(0) = t0;
    x(1) = t1;
    return x;
}

template <typename T>
MultiVariatePoint<T> MultiVariatePoint<T>::toMonoVariatePoint(T t)
{
    return MultiVariatePoint<T>(1,t);
}

template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator=(const MultiVariatePoint<T> & nu)
{
    if (this != &nu)
    {
        m_d = nu.getD();
        if (m_d != 0)
        {
            delete[] m_nu;
            m_nu = new T[m_d];
            m_waitingTime = 0;
            m_alpha = 0;
            memcpy(m_nu,nu.m_nu,sizeof(T)*nu.getD());
        }
    }
    return *this;
}

template <typename T>
MultiVariatePoint<T>& MultiVariatePoint<T>::operator+=(const MultiVariatePoint<T> & nu)
{
    for (int i=0; i<nu.getD(); i++)
        (*this)(i) += nu(i);
    return *this;
}

template <typename T>
std::ostream & operator<<(std::ostream &out, const MultiVariatePoint<T> &nu)
{
    out << "(";
    for (int k=0; k<int(nu.getD()-1); k++)
        cout << nu(k) << ",";
    out << nu(nu.getD()-1) << ")";
    return out;
}


template <typename T>
bool operator<(MultiVariatePoint<T> const &nu1, MultiVariatePoint<T> const &nu2)
{
    int k=0;
    while (nu1(k)==nu2(k) && k<int(nu1.getD()-1))
        k++;
    return (nu1(k)<nu2(k));
}

template <typename T>
bool operator==(MultiVariatePoint<T> const &nu1 , MultiVariatePoint<T> const &nu2)
{
    if (nu1.getD() != nu2.getD()) return false ;
    int i = 0;
    while (i < nu1.getD())
    {
        if (nu1(i) != nu2(i)) return false ;
        i++ ;
    }
    return true ;
}

#endif
