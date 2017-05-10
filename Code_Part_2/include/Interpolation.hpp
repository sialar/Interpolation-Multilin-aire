#ifndef INTERPOLATION
#define INTERPOLATION

#include <iostream>
#include <vector>
#include <list>
#include <cmath>

#include "BinaryTree.hpp"

using namespace std;

class Interpolation
{
    private:
        vector<double> m_points; // ordre: -1,1,0,-0.5,0.5,-0.75,-0.25,0.25,0.75,....
        vector<double> m_alphaTab; // m_alphaTab[i] correspond a alpha de ieme point dans m_points
        vector<int> m_path; // m_path[i] correspond au ieme élément pris dans l'ordre imposé par l'algo AI
                            // m_points[m_path[i]] la valeur du ieme point dans path
        list<int> m_nextPoints; // contient les indices des prochains points potentiels
        BinaryTree* m_tree;

    public:

        Interpolation(int size);
        ~Interpolation();

        double piecewiseFunction_1D(int k, double y);
        double computeLastAlphaI(int i);
        void updateNextPoints(int i);
        bool indiceInPath(int index);
        int getIndice(double l);

        double interpolation_iterative(double y, int k, bool debug);
        static double g(double y) { return y*y; };

        void displayPath();
        void displayAlphaTab();
        void displayNextPoints();
};
#endif
