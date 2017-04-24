#include "../include/IndiceND.hpp"

IndiceND::IndiceND()
{
    m_d = 0;
}

IndiceND::IndiceND(int d, int val)
{
    m_d = d;
    if (m_d != 0)
    {
        m_nu = new int[m_d];
        for (int i=0; i<m_d; i++)
            m_nu[i] = val;
    }
}

IndiceND::IndiceND(const IndiceND& nu)
{
    *this=nu;
}

IndiceND::~IndiceND()
{
    if (m_nu != NULL)
        delete[] m_nu;
}

IndiceND& IndiceND::operator=(const IndiceND & nu)
{
    if (this != &nu)
    {
        m_d = nu.getD();
        if (m_d != 0)
        {
            delete[] m_nu;
            m_nu = new int[m_d];
            memcpy(m_nu,nu.m_nu,sizeof(int)*nu.getD());
        }
    }
    return *this;
}

IndiceND& IndiceND::operator+=(const IndiceND & nu)
{
    for (int i=0; i<nu.getD(); i++)
        (*this)(i) += nu(i);
    return *this;
}

bool operator<(IndiceND const &nu1, IndiceND const &nu2)
{
    for (int k=0; k<nu1.getD(); k++)
        if (nu1(k)>nu2(k));
            return false;
    return true;
}

std::ostream & operator<<(std::ostream &out, const IndiceND &nu)
{
    out << "(";
    for (int k=0; k<int(nu.getD()-1); k++)
        cout << nu(k) << ",";
    out << nu(nu.getD()-1) << ")";
    return out;
}

bool operator==(IndiceND const &nu1 , IndiceND const &nu2)
{
    if (nu1.getD() != nu2.getD())
        return false ;
    else
    {
        bool equals = true ;
        int i = 0;
        while (equals && (i < nu1.getD()))
        {
            if (nu1(i) != nu2(i))
                equals = false ;
            i++ ;
        }
        return equals ;
    }
}
