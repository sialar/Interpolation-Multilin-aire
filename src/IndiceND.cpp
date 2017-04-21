#include "../include/IndiceND.hpp"

IndiceND::IndiceND(int d, int val)
{
    m_d = d;
    m_i.resize(d,val);
}

IndiceND::IndiceND(const IndiceND& nu)
{
    copy(nu);
}

IndiceND::~IndiceND()
{
}

void IndiceND::display()
{
    cout << "(";
    for (int k=0; k<int(m_i.size()-1); k++)
        cout << m_i[k] << ",";
    cout << m_i.back() << ")" << endl;
}

void IndiceND::copy(const IndiceND& nu)
{
    m_d = nu.m_d;
    m_i.resize(m_d);
    for (int k=0; k<m_d; k++)
        m_i[k] = nu.m_i[k];
}

int IndiceND::operator()(int d) const {
    if ((d >= 0) && (d <= m_d-1))
        return m_i[d];
    else
        cerr << "Erreur: element indisponible" << endl;
    return -1;
}

IndiceND& IndiceND::operator =(const IndiceND & nu){
    m_d = nu.getD();
    m_i.resize(m_d);
    for (int i=0; i<m_d; i++)
        m_i[i] = nu(i);
    return *this;
}

std::ostream & operator<<(std::ostream &out, const IndiceND &nu)
{
    out << "(";
    for (int k=0; k<int(nu.getD()-1); k++)
        cout << nu(k) << ",";
    out << nu(nu.getD()-1) << ")" << endl;
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
