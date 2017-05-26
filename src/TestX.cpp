#include <iostream>
#include <string>
#include <time.h>
#include "../include/Utils.hpp"
#include "../include/MultiVariatePoint.hpp"

using namespace std;

void compareAndDisplay(MultiVariatePoint<string> s1, MultiVariatePoint<string> s2)
{
    cout << s1;
    if (Utils::equals(s1,s2)) cout << " == ";
    else cout << " != ";
    cout << s2 << endl;
}
int main( int argc, char* argv[] )
{
    srand (time(NULL));

    MultiVariatePoint<string> s1(2,"010");
    MultiVariatePoint<string> s2(2,"010");
    MultiVariatePoint<string> s3(2,"011");

    compareAndDisplay(s1,s2);
    compareAndDisplay(s1,s3);
    compareAndDisplay(s2,s3);
    Utils::storeDichotomySequenceInFile(33);
    return 0;
}
