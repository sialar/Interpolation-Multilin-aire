#include <iostream>
#include <string>
#include <array>
#include <time.h>
#include <limits>
#include <stdio.h>
#include "../include/Tucker/TuckerApproximation.hpp"
#include "../include/Utils.hpp"

using namespace std;


#define BUFSIZE 128

int parse_output(void) {
    string rep = Utils::projectPath + "Tucker/LIVRAISON_THESE_Paris6_MOX_test";
    string s = "cd " + rep + "\npython eval.py macro_fission0 macro_scattering000 0 0 0 0 0";
    char *cmd = new char[s.length() + 1];
    strcpy(cmd, s.c_str());

    char buf[BUFSIZE];
    FILE *fp;

    if ((fp = popen(cmd, "r")) == NULL) {
        printf("Error opening pipe!\n");
        return -1;
    }

    while (fgets(buf, BUFSIZE, fp) != NULL) {
        // Do whatever you want here...
        printf("OUTPUT: %s", buf);
    }

    if(pclose(fp))  {
        printf("Command not found or exited with error status\n");
        return -1;
    }

    return 0;
}

int main( int argc, char* argv[] )
{
    int a = parse_output();
    return 0;
}
