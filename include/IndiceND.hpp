#ifndef INDICEND
#define INDICEND

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>

using namespace std;

class IndiceND
{
    private:
        int m_d;
        vector<int> m_i;
        bool m_inPath;

    public:
        IndiceND(int d, int val=0);
        IndiceND(const IndiceND& nu);
        ~IndiceND();

        int getD() const { return m_d; };
        vector<int> getValues() { return m_i; };
        void setValue(int k, int val) { m_i[k] = val; };

        void copy(const IndiceND& nu);

        bool isInPath() { return m_inPath; };
        void setInPath(bool val) { m_inPath = val; };

        void display();

        int operator()(int d) const;
        IndiceND& operator=(IndiceND const &nu);
};

std::ostream & operator<<(std::ostream &out, const IndiceND &nu) ;
bool operator==(IndiceND const &nu1 , IndiceND const &nu2);

#endif
