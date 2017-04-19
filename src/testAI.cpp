#include <iostream>
#include <string>
#include <time.h>
#include "../include/LagrangeInterpolation2D.hpp"
#include "../include/Utils.hpp"

using namespace std;

void displayTestPoints(vector<double> vx, vector<double> vy);
vector<double> displayGRealValues(LagrangeInterpolation2D* interp, vector<double> vx, vector<double> vy);


int main( int argc, char* argv[] )
{
    srand (time(NULL));

    // Get the size of the data points sequence from argument
    int n = (argc > 1) ? stoi(argv[1]) : Utils::randomValue(10,100);
    int m = (argc > 2) ? stoi(argv[2]) : Utils::randomValue(10,100);

    int nbTestPointsX = (argc > 3) ? stoi(argv[3]) : Utils::randomValue(2,5);
    int nbTestPointsY = (argc > 4) ? stoi(argv[4]) : Utils::randomValue(2,5);
    LagrangeInterpolation2D* interp = new LagrangeInterpolation2D(n, m, -1);

    // Evaluation de la fonction g (points choisis aléatoirement entre -1 et 1)
    vector<double> testPointsX, testPointsY, realValues, estimate;
    Utils::separateur();

    testPointsX.resize(nbTestPointsX);
    for (int i=0; i<nbTestPointsX; i++)
        testPointsX[i] = Utils::randomValue(-1,1);
    testPointsY.resize(nbTestPointsY);
    for (int i=0; i<nbTestPointsY; i++)
        testPointsY[i] = Utils::randomValue(-1,1);
    displayTestPoints(testPointsX,testPointsY);
    Utils::separateur();
    realValues = displayGRealValues(interp,testPointsX,testPointsY);
    Utils::separateur();



    interp->setPointsX(Utils::createLejaSequence(n));
    interp->setPointsY(Utils::createLejaSequence(m));
    interp->showAlphaTab();

    interp->buildPathWithAIAlgo(n,m,true);

    cout << "   - Chemin obtenu avec l'algo AI:" << endl;
    interp->showPath();
    Utils::separateur();


    // Test de l'interpolation en utilisant la sequence de Leja
    double val = 0;
    cout << endl << "   - Calcul par interpolation bilinéaire: (evaluation de ĝ en " << nbTestPointsX*nbTestPointsY << " points de test)" << endl;
    interp->computeAllAlphaNu();
    for (int i=0; i<int(testPointsX.size()); i++)
    {
        for (int j=0; j<int(testPointsY.size()); j++)
        {
            val = interp->lagrangeInterpolation_2D_iterative(testPointsX[i],testPointsY[j]);
            estimate.push_back(val);
            cout << val << " ";
        }
        cout << endl;
    }
    cout << "   - Erreur quadratique moyenne: " << Utils::squareError(realValues,estimate) << endl;
    Utils::separateur();

    return 0;
}


void displayTestPoints(vector<double> vx, vector<double> vy)
{
  cout << "   - Sequence des points de test: (" << vx.size() << " suivant l'axe des x et " <<
  vy.size() << " suivant l'axe des y):" << endl;
  cout << "         + Suivant la direction x: ";
  for (int i=0; i<int(vx.size()); i++)
  {
    vx[i] = Utils::randomValue(-1,1);
    cout << vx[i] << " ";
  }
  cout << endl << "         + Suivant la direction y: ";
  for (int j=0; j<int(vy.size()); j++)
  {
    vy[j] = Utils::randomValue(-1,1);
    cout << vy[j] << " ";
  }
  cout << endl ;
}
vector<double> displayGRealValues(LagrangeInterpolation2D* interp, vector<double> vx, vector<double> vy)
{
  vector<double> realValues;
  double val = 0;
  cout << "   - Calcul direct (evaluation de g en " << vx.size()*vy.size() << " points):" << endl;
  for (int i=0; i<int(vx.size()); i++)
  {
    for (int j=0; j<int(vy.size()); j++)
    {
      val = interp->g(vx[i],vy[j]);
      realValues.push_back(val);
      cout << val << " ";
    }
    cout << endl;
  }
  return realValues;
}
