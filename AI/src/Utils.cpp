#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;
double Utils::m_precision = numeric_limits<double>::digits10+1;


/******************************************************************************/
/************************ Fonctions d'affichage *******************************/
void Utils::displayValues(vector<double> values)
{
    cout << "[ " << m_precision << values[0];
    for (size_t i=1; i<values.size(); ++i)
        cout << " ; " <<  m_precision << values[i];
    cout << " ]" << endl;
}

/******************************************************************************/
/****************** Opérations sur les vecteur de double **********************/
double Utils::min_elt(vector<double> v)
{
    // Retourner l'élément max de v
    if (v.empty()) return 0;
    double m = v[0];
    for (double x : v)
        if (x<m)
            m = x;
    return m;
}
double Utils::max_elt(vector<double> v)
{
    // Retourner l'élément min de v
    if (v.empty()) return 0;
    double m = v[0];
    for (double x : v)
        if (x>m)
            m = x;
    return m;
}
bool absLess(double x, double y)
{
    // Comparer la valeur absolue de deux valeurs
    return abs(x) < abs(y);
}
double Utils::norm(vector<double> x, int p)
{
    // Calculer la norme d'un vecteur x
    // si p=0 --> norme infinie
    // sinon --> norme 2
    if (!p) return abs(*max_element(x.begin(), x.end(), absLess));
    else
    {
        double res = 0.0;
        for (int i=0; i<int(x.size()); i++)
            res += pow(x[i],2);
        return sqrt(res);
    }
}
vector<double> Utils::diff(vector<double> x, vector<double> y)
{
    // Calculer la différence entre 2 vecteurs x et y
    vector<double> res;
    for (int i=0; i<int(x.size()); i++)
        res.push_back(x[i]-y[i]);
    return res;
}

/******************************************************************************/
/**************** Opérations de changement de variables ***********************/
double Utils::convertToDefaultDomain(double a, double b, double x)
{
    // Changement de variable [a,b] --> [-1,1]
    return (2*x)/(b-a) + 1 - (2*b)/(b-a);
}
double Utils::convertToFunctionDomain(double a, double b, double x)
{
    // Changement de variable [-1,1] --> [a,b]
    return x*(b-a)/2 + (a+b)/2;
}

/******************************************************************************/
/**************** Opérations sur des chaines de caractères ********************/
string Utils::replace(string strs, string str_old, string str_new)
{
  // Remplacer, dans la chaine de caractères strs, les sous-chaines str_old par str_new
  // Recharche de str_old
  size_t found = strs.find(str_old);
  while (found!=string::npos)
  {
    // Tant qu'on trouve str_old on poursuis le remplacement
    strs.replace(found, str_old.length(), str_new);
    // Recharche de str_old
    found = strs.find(str_old);
  }
  // Retourner la nouvelle chaine
  return strs;
}
string Utils::eraseExtraSpaces(string strs)
{
    // Supprimer les espaces de séparation inutiles dans strs
    // En sortie strs est une chaine de caractères contenant des réels séparés par un seul espace
    string s = replace(strs, "\t", " ");
    return replace(s, "  ", " ");
}
bool Utils::strInVector(string required, vector<string> vec)
{
    // Vérifier si la chaine de caractère required est un élément du vecteur vec
    for (string s : vec)
        if (required.compare(s)==0)
            return true;
    return false;
}
string Utils::vector2str(vector<double> x)
{
    // Convertir un vecteur de double en chaine de caractères
    string s = "(" + to_string(x[0]);
    for (int i=1; i<int(x.size()); i++)
        s = s + ", " + to_string(x[i]);
    s = s + ") ";
    return s;
}
vector<double> Utils::str2vector(string line)
{
    // Convertir une chaine de caractères en un vecteur de double
    // On commence par supprimer les espace inutiles
    string str = Utils::eraseExtraSpaces(line);
    stringstream ss(str);
    vector<double> data;
    string word;
    // On parcours la chaine mot par mot
    while (getline(ss, word, ' '))
        // Conversion du mot courant en double et ajout au vecteur data
        data.push_back(stod(word));
    return data;
}

/******************************************************************************/
/******** Foncions utiles pour la création des points d'interpolation *********/
double Utils::randomValue(double a, double b)
{
  return (rand()/(double)RAND_MAX) * (b-a) + a;
}
MultiVariatePoint<double> Utils::createRandomMultiVariatePoint(int d)
{
    MultiVariatePoint<double> point(d,0,0);
    for (int i=0; i<d; i++)
        point(i) = Utils::randomValue(-1,1);
    return point;
}
vector<double> Utils::createUniformSequence(int nbPoints)
{
    // Construire une séquence uniforme de nbPoints points
    // Voir page 21 de l'article High-dimensional adaptive sparse polynomial interpolation
    // and applications to parametric PDEs
    double sum = 0;
    vector<double> points;
    vector<double> binary_exp;
    int stop = (nbPoints%2) ? (nbPoints-1)/2 : (nbPoints-1)/2+1;
    points.resize(nbPoints);

    // 3 premier points de la séquence uniforme
    points[0] =  1;
    points[1] = -1;
    points[2] =  0;
    for (int k=1; k<stop; ++k)
    {
        sum = 0;
        binaryExpansion(k,binary_exp);
        for (size_t j=0; j<binary_exp.size(); ++j)
            sum += binary_exp[j] / pow(2,j);

        points[2*k+1] = 0.5 * sum;
        points[2*k+2] = - points[2*k+1];
    }
    return points;
}
vector<double> Utils::createChebychevSequence(int nbPoints)
{
    // Créer une séquence de nbPoints points de Chebychev
    vector<double> points;
    points.resize(nbPoints);
    for (int i=0; i<nbPoints; i++)
        points[i] = cos(i*M_PI/(nbPoints-1));
    return points;
}
vector<double> Utils::createLejaSequence(int nbPoints)
{
    // Créer une séquence de nbPoints points de Leja
    // Les points sont choisis parmi la grille fine m_1dGrid
    m_1dGrid = createChebychevSequence(10001);
    vector<double> points;
    // Premier point
    points.push_back(1);
    double newPoint;
    for (int i=1; i<nbPoints; i++)
    {
        // Tant qu'on ne dépasse pas la taille souhaitée
        while (int(points.size()) < nbPoints)
        {
            // Calcul du nouveau point de Leja
            newPoint = computeNewLejaPointFromSequence(points);
            points.push_back(newPoint);
        }
    }
    return points;
}
void Utils::storeLejaSequenceInFile(int nbPoints)
{
  // Enregister la séquence de Leja de longeur nbPoints dans le fichier data/leja_sequence.dat
  ofstream file("AI/data/leja_sequence.dat", ios::out | ios::trunc);
  if(file)
  {
    vector<double> lejaSeq = createLejaSequence(nbPoints);
    for (int i=0; i<nbPoints; ++i)
    file << lejaSeq[i] << endl;
    file.close();
  }
  else
  cerr << "Error while opening the file!" << endl;
}
vector<double> Utils::loadLejaSequenceFromFile(int nbPoints)
{
  // Charger la séquence de Leja de longeur length depuis le fichier data/leja_sequence.dat
  // Pour éviter de recalculer les points de Leja à chaque éxecution
  ifstream file("AI/data/leja_sequence.dat", ios::in);
  vector<double> lejaSeq;
  string line;
  int line_index = 0;
  if(file)
  {
    while (getline(file,line) && line_index<nbPoints)
    {
      lejaSeq.push_back(stof(line));
      line_index++;
    }
    file.close();
  }
  else
  cerr << "Error while opening the file! " << endl;
  return lejaSeq;
}
void Utils::binaryExpansion(int number, vector<double>& binary_exp)
{
  // Calculer l'expansion binaire d'un entier number
  // Utile pour la création de sequence uniforme
  // number = \sum_{j=0}^n epsilon_j 2^j
  // en sortie binary_exp = {epsilon_j}_j
  binary_exp.clear();
  int temp = number;
  while (temp > 0)
  {
    binary_exp.push_back(temp % 2);
    temp = temp / 2;
  }
}
bool Utils::isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold)
{
    // Renvoie True si la distance de y à un point de seq est inférieure à threshold
    for (int i=0; i<int(seq.size()); i++)
        if (abs(y-seq[i]) <= threshold)
            return true;
    return false;
}

double Utils::computeNewLejaPointFromSequence(vector<double> seq)
{
    // Calculer le nouveau point de Leja dans la liste de points seq
    int argmax = -1;
    double prod = 1;
    double max = numeric_limits<double>::min();
    for (int k=0; k<int(m_1dGrid.size()); k++)
    {
        prod = 1;
        // si le point courant est très proche d'un point de la séquence de Leja courante,
        // on le prend pas en compte
        // sinon on verifie s'il est le plus loin des points de seq
        if (!isTooCloseToOneLejaPoint(m_1dGrid[k],seq,1e-10))
        {
            for (int i=0; i<int(seq.size()); i++)
                prod *= abs(m_1dGrid[k] - seq[i]);
            if (prod > max)
            {
                max = prod;
                argmax = k;
            }
        }
    }
    return m_1dGrid[argmax];
}

/******************************************************************************/
/************** Foncions de comparaison de points multivariés *****************/
bool Utils::equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2)
{
    // Comparer 2 points multivariés
    if (nu1.getD() != nu2.getD()) return false;
    for (int i=0; i<nu1.getD(); i++)
        if (nu1(i).compare(nu2(i))!=0)
            return false;
    return true;
}

/******************************************************************************/
/********** Foncions utils pour la lecture des entrées du programme ***********/
int Utils::chooseDimensionD(int argc, char* argv[], int argNum)
{
    // Lecture depuis les arguments d'execution
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    // Lecture interactive
    while (dim < 0)
    {
        cout << " - Choose the space dimension d: ";
        cin >> dim;
    }
    return dim;
}
int Utils::chooseDimensionN(int argc, char* argv[], int argNum)
{
    // Lecture depuis les arguments d'execution
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    // Lecture interactive
    while (dim < 0)
    {
        cout << " - Choose the space dimension n: ";
        cin >> dim;
    }
    return dim;
}
int Utils::chooseNbTestPoints(int argc, char* argv[], int argNum)
{
  // Lecture depuis les arguments d'execution
  if (argc > argNum) return stoi(argv[argNum]);
  int nbTestPoints = -1;
  // Lecture interactive
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}
int Utils::chooseMaxIteration(int argc, char* argv[], int argNum)
{
  // Lecture depuis les arguments d'execution
  if (argc > argNum) return stoi(argv[argNum]);
  int maxIteration = -1;
  // Lecture interactive
  while (maxIteration < 0)
  {
    cout << " - Choose the maximum number of iteration : ";
    cin >> maxIteration;
  }
  return maxIteration;
}
MultiVariatePoint<int> Utils::chooseMethods(int dim)
{
    MultiVariatePoint<int> methods(dim,0,-1);
    for (int i=0; i<dim; i++)
    {
        while (methods(i)!=0 && methods(i)!=1 && methods(i)!=2)
        {
            cout << " - Choose the method of interpolation in direction [" << i << "]: " << endl;
            cout << "\t - 0: Using lagrange polynomial functions and leja points: " << endl;
            cout << "\t - 1: Using piecewise functions and middle points: " << endl;
            cout << "\t - 2: Using quadratic functions and middle points: " << endl << " - ";
            cin >> methods(i);
        }
    }
    return methods;
}
