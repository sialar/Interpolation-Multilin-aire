#include <iostream>
#include <string>
#include "../include/BinaryTree.hpp"
#include "../include/Utils.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
    int depth = (argc>1) ? stoi(argv[1]) : Utils::randomValue(0,5);
    double target = (argc>2) ? stof(argv[2]) : Utils::randomValue(0,1);
    BinaryTree* tree = new BinaryTree();
    tree->initTree(depth);
    tree->displayBinaryTree();

    double res_sup, res_inf;
    tree->searchNode(target,&res_sup,&res_inf);
    cout << target << " found: (" << res_inf << "," << res_sup << ")" << endl;
    return 0;
}
