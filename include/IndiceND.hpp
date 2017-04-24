#ifndef INDICEND
#define INDICEND

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstring>

using namespace std;

class IndiceND
{
    private:
        int m_d;
        int* m_nu = NULL;
        bool m_inPath;

    public:
        IndiceND();
        IndiceND(int d, int val=0);
        IndiceND(const IndiceND& nu);
        ~IndiceND();

        int getD() const { return m_d; };

        bool isInPath() { return m_inPath; };
        void setInPath(bool val) { m_inPath = val; };

        int &operator()(int d) const {return m_nu[d]; };
        IndiceND& operator+=(const IndiceND & v);
        IndiceND& operator=(IndiceND const &nu);
};
bool operator< (IndiceND const &nu1, IndiceND const &nu2);
std::ostream & operator<<(std::ostream &out, const IndiceND &nu) ;
bool operator==(IndiceND const &nu1 , IndiceND const &nu2);

#endif
