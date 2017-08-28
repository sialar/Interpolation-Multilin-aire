#include "../include/Utils.hpp"

vector<double> Utils::m_1dGrid;
double Utils::m_precision = numeric_limits<double>::digits10+1;

void Utils::separateur()
{
    for (int i=0; i<206; i++) cout << "=";
    cout << endl;
}

// Afficher les éléments d'un vecteur
void Utils::displayValues(vector<double> values)
{
    cout << "[ " << m_precision << values[0];
    for (size_t i=1; i<values.size(); ++i)
        cout << " ; " <<  m_precision << values[i];
    cout << " ]" << endl;
}

// Changement de variable [a,b] --> [-1,1]
double Utils::convertToDefaultDomain(double a, double b, double x)
{
    return (2*x)/(b-a) + 1 - (2*b)/(b-a);
}

// Changement de variable [-1,1] --> [a,b]
double Utils::convertToFunctionDomain(double a, double b, double x)
{
    return x*(b-a)/2 + (a+b)/2;
}

// Retourner une valeur entre a et b aléatoirement
double Utils::randomValue(double a, double b)
{
  return (rand()/(double)RAND_MAX) * (b-a) + a;
}

// créer un vecteur aléatoirement
MultiVariatePoint<double> Utils::createRandomMultiVariatePoint(int d)
{
    MultiVariatePoint<double> point(d,0,0);
    for (int i=0; i<d; i++)
        point(i) = Utils::randomValue(-1,1);
    return point;
}

// Comparer la valeur absolue de deux valeurs
bool absLess(double x, double y)
{
    return abs(x) < abs(y);
}

// Calculer la norme d'un vecteur x
// si p=0 --> norme infinie
// sinon --> norme 2
double Utils::norm(vector<double> x, int p)
{
    if (!p) return abs(*max_element(x.begin(), x.end(), absLess));
    else
    {
        double res = 0.0;
        for (int i=0; i<int(x.size()); i++)
            res += pow(x[i],2);
        return sqrt(res);
    }
}

// Calculer la différence entre 2 vecteurs
vector<double> Utils::diff(vector<double> x,vector<double> y)
{
    vector<double> res;
    for (int i=0; i<int(x.size()); i++)
        res.push_back(x[i]-y[i]);
    return res;
}

// Retourner le plus grand coeficient en valeur absolue
double Utils::maxAbsValue(vector<double> vec)
{
    vector<double> vec_copy;
    for (int i=0; i<int(vec.size()); i++)
        vec_copy.push_back(abs(vec[i]));
    return max_elt(vec_copy);
}

// Renvoie True si la distance de y à un point de seq est inférieure à threshold
bool Utils::isTooCloseToOneLejaPoint(double y, vector<double> seq, double threshold)
{
    for (int i=0; i<int(seq.size()); i++)
        if (abs(y-seq[i]) <= threshold)
            return true;
    return false;
}

// Calculer le nouveau point de Leja dans la liste de points seq
double Utils::computeNewLejaPointFromSequence(vector<double> seq)
{
    int argmax = -1;
    double prod = 1;
    double max = numeric_limits<double>::min();
    for (int k=0; k<int(m_1dGrid.size()); k++)
    {
        prod = 1;
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

// Décomposer le le nombre number en binaire
void Utils::binaryDecomposition(int number, vector<double>& binary_decomp)
{
  binary_decomp.clear();
  int temp = number;
  while (temp > 0)
  {
    binary_decomp.push_back(temp % 2);
    temp = temp / 2;
  }
}

// Enregister la séquence de Leja de longeur length dans le fichier data/leja_sequence.dat
void Utils::storeLejaSequenceInFile(int length)
{
    ofstream file("AI/data/leja_sequence.dat", ios::out | ios::trunc);
    if(file)
    {
        vector<double> lejaSeq = createLejaSequence(length);
        for (int i=0; i<length; ++i)
              file << lejaSeq[i] << endl;
        file.close();
    }
    else
        cerr << "Error while opening the file!" << endl;
}

// Charger la séquence de Leja de longeur length depuis le fichier data/leja_sequence.dat
// Pour éviter de recalculer les points de Leja à chaque éxecution
vector<double> Utils::loadLejaSequenceFromFile(int length)
{
    ifstream file("AI/data/leja_sequence.dat", ios::in);
    vector<double> lejaSeq;
    string line;
    int line_index = 0;
    if(file)
    {
        while (getline(file,line) && line_index<length)
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

// Vérifier si la chaine de caractère required est un élément du vecteur vec
bool Utils::strInVector(string required, vector<string> vec)
{
    for (string s : vec)
        if (required.compare(s)==0)
            return true;
    return false;
}

// Convertir un vecteur de double en chaine de caractères
string Utils::vector2str(vector<double> x)
{
    string s = "(" + to_string(x[0]);
    for (int i=1; i<int(x.size()); i++)
        s = s + ", " + to_string(x[i]);
    s = s + ") ";
    return s;
}

// Convertir une chaine de caractères en un vecteur de double
vector<double> Utils::str2vector(string line)
{
    string str = Utils::eraseExtraSpaces(line);
    stringstream ss(str);
    vector<double> data;
    string word;
    while (getline(ss, word, ' '))
        data.push_back(stod(word));
    return data;
}

// Calculer l'erreur relative absolue
double Utils::relativeError(vector<double> f,vector<double> f_tilde)
{
    int n = min(int(f.size()),int(f_tilde.size()));
	vector<double> e;
    for (int k=0; k<n; k++)
        e.push_back((f_tilde[k] - f[k]) * pow(10,5) / maxAbsValue(f));
    return maxAbsValue(e);
}

// Calculer l'erreur quadratique moyenne
double Utils::computeMseError(vector<double> f,vector<double> f_tilde)
{
    double e = 0.0;
    int n = min(int(f.size()),int(f_tilde.size()));
    for (int k=0; k<n; k++)
        e += pow((f_tilde[k] - f[k]) * pow(10,5),2);
    return sqrt(e/n);
}

// Construire une séquence uniforme de nbPoints points
vector<double> Utils::createUniformSequence(int nbPoints)
{
    double sum = 0;
    vector<double> points;
    vector<double> binary_decomp;
    int stop = (nbPoints%2) ? (nbPoints-1)/2 : (nbPoints-1)/2+1;
    points.resize(nbPoints);

    points[0] =  1;
    points[1] = -1;
    points[2] =  0;
    for (int k=1; k<stop; ++k)
    {
        sum = 0;
        binaryDecomposition(k,binary_decomp);
        for (size_t j=0; j<binary_decomp.size(); ++j)
            sum += binary_decomp[j] / pow(2,j);

        points[2*k+1] = 0.5 * sum;
        points[2*k+2] = - points[2*k+1];
    }
    return points;
}

// Créer une séquence de nbPoints points de Chebychev
vector<double> Utils::createChebychevSequence(int nbPoints)
{
    vector<double> points;
    points.resize(nbPoints);
    for (int i=0; i<nbPoints; i++)
        points[i] = cos(i*M_PI/(nbPoints-1));
    return points;
}

// Créer une séquence de nbPoints points de Leja
vector<double> Utils::createLejaSequence(int nbPoints)
{
    m_1dGrid = createChebychevSequence(10001);
    vector<double> points;
    points.push_back(1);
    double newPoint;
    for (int i=1; i<nbPoints; i++)
    {
        while (int(points.size()) < nbPoints)
        {
            newPoint = computeNewLejaPointFromSequence(points);
            points.push_back(newPoint);
        }
    }
    return points;
}

// Comparer 2 points multivariés
bool Utils::equals(MultiVariatePoint<string> nu1, MultiVariatePoint<string> nu2)
{
    if (nu1.getD() != nu2.getD()) return false;
    for (int i=0; i<nu1.getD(); i++)
        if (nu1(i).compare(nu2(i))!=0)
            return false;
    return true;
}

// Remplacer, dans la chaine de caractères strs, la sous-chaine str_old par str_new
string Utils::replace(string strs, string str_old, string str_new)
{
  size_t found = strs.find(str_old);
  while (found!=string::npos)
  {
    strs.replace(found, str_old.length(), str_new);
    found = strs.find(str_old);
  }
  return strs;
}

// Supprimer les espaces de séparation inutiles dans strs
// En sortie strs est une chaine de caractères contenant des réels séparés par un seul espace
string Utils::eraseExtraSpaces(string strs)
{
    string s = replace(strs, "\t", " ");
    return replace(s, "  ", " ");
}

// Retourner l'élément max de v
double Utils::min_elt(vector<double> v)
{
    if (v.empty()) return -11;
    double m = v[0];
    for (double x : v)
        if (x<m)
            m = x;
    return m;
}

// Retourner l'élément min de v
double Utils::max_elt(vector<double> v)
{
    if (v.empty()) return -11;
    double m = v[0];
    for (double x : v)
        if (x>m)
            m = x;
    return m;
}

int Utils::chooseDimensionD(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension d: ";
        cin >> dim;
    }
    return dim;
}

int Utils::chooseDimensionN(int argc, char* argv[], int argNum)
{
    if (argc > argNum) return stoi(argv[argNum]);
    int dim = -1;
    while (dim < 0)
    {
        cout << " - Choose the space dimension n: ";
        cin >> dim;
    }
    return dim;
}

int Utils::chooseNbTestPoints(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int nbTestPoints = -1;
  while (nbTestPoints < 0)
  {
    cout << " - Choose the number ot test points : ";
    cin >> nbTestPoints;
  }
  return nbTestPoints;
}

int Utils::chooseMaxIteration(int argc, char* argv[], int argNum)
{
  if (argc > argNum) return stoi(argv[argNum]);
  int maxIteration = -1;
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
