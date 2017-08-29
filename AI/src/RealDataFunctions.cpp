#include "../include/RealDataFunctions.hpp"

RealDataFunctions::RealDataFunctions(int d, int n, string fileName) : Functions(d,n)
{
	m_n = n;
	m_d = d;
	m_fileName = fileName;
	readDataFromFile();
	// Dans cet exemple, on va approcher la section efficace associé à la réaction totale de type MOX
	// Ceci va servir dans la méthode evaluate, pour donner la valeur supposée "exacte" en n'importe quel point x
	// N'est utile que pour l'approximation des sections efficaces quand on suppose que la méthode de Tucker
	// nous donne les valeurs "exactes"
	m_tuckerApprox = make_shared<TuckerApproximation>("MOX","macro_totale0");
}

vector<double> RealDataFunctions::evaluate(MultiVariatePoint<double> x)
{
		// En entrée x appartient à [-1,1]
		// On se met dans l'éspace réel de f par changement de variable, avant de l'évaluer en un point x.
		// Dans chaque direction: x : [-1,1] --> [a,b]_j domaine réel de f sur la direction j
		for (int i=0; i<m_d; i++)
				x(i) = Utils::convertToFunctionDomain(m_parametersDomain[i][0], m_parametersDomain[i][1], x(i));

		// Pour respecter le type de retour dans le cas général, on retourne un vecteur
		// même quand f est à valeurs réelles
		vector<double> res(1,m_tuckerApprox->evaluate(x,"macro_totale0"));
    return res;
}

void RealDataFunctions::readDataFromFile()
{
		// Ouverture du fichier "AI/data/"m_fileName
    ifstream file("AI/data/" + m_fileName, ios::in);
    if(file)
  	{
				// Contient la ligne courante du fichier
	    	string line;
				// Ligne courante sous forme de vecteur de double
	      vector<double> data;
				// max et min permettent de construire le vrai domaine de f (non pas [-1,1] qui est le domaine par défaut
				// utilisé dans l'algorithme AI). Une fois la lecture du fichier terminée, on sauvegarde la valeur max et min
				// sur chaque direction
	      vector<double> max(m_d,-numeric_limits<double>::max());
	      vector<double> min(m_d,numeric_limits<double>::max());
				// p = (p1, .., pd) : point de référence
	      MultiVariatePoint<double> p(m_d,0,0);

				// Parcourir le fichier ligne par ligne
				// on note li la ligne i: li = (indice, p1, .., pd, valeur_exacte)
	      while (getline(file, line))
	      {
						// Conversion de chaine de caractère line en vecteur de double
	          data = Utils::str2vector(line);

						// Pour vérifier que les dimensions données en entrée sont conformes aux dimensions présente dans
						// le fichier de données. Si ce n'est pas le cas, arrêt de l'éxecution
						if (int(data.size()) != m_d + m_n + 1)
						{
								cerr << "Choix des dimensions (n et d) et contenu du fichier incompatibles" << endl;
								cout << data.size() << " " << m_d + m_n + 1 << endl;
								exit(1);
						}

						// Lecture du point de référence (de taille m_d)
        		for (int i=0; i<m_d; i++)
	              p(i) = data[i+1];
						// Ajout du point à la liste des points de référence
        		m_referencePts.push_back(p);
						// Mise à jour de max et min pour le calcul du domaine réel de f
			      for (int i=0; i<m_d; i++)
			      {
			      		if (p(i) > max[i]) max[i] = p(i);
			          if (p(i) < min[i]) min[i] = p(i);
			      }
						// Lecture de la valeur de f au point de référence courant p
						// On prend les derniers éléments dans le vecteur data (on enleve indice et les pi)
						// On général, on une seule valeur sauf pour le cas des fonctions à valeurs vectorielles
						for (int i=0; i<m_d+1; i++)
								data.erase(data.begin());
						// Ajout de la valeur à la liste des valeurs de référence de f
						m_exactValues.push_back(data);
	      }
				// Construction des domaines réels de f, une fois la lecture terminée
	      for (int i=0; i<m_d; i++)
	      {
						m_parametersDomain[i][0] = min[i];
						m_parametersDomain[i][1] = max[i];
	      }
				// Conversion des point de référence par changement de variable direction par direction [a,b]_j --> [-1,1]
				// Pourquoi ?
				// --> L'algorithme d'interpolation adaptative est implémenté pour des fonctions définies sur [-1,1]^d
				// Pour prendre en compte n'importe quel ensemble de départ, on procède par changement de variable.
	      for (int i=0; i<int(m_referencePts.size()); i++)
	          for (int j=0; j<m_d; j++)
	              m_referencePts[i](j) = Utils::convertToDefaultDomain(min[j], max[j], m_referencePts[i](j));
	      file.close();
	  }
	  else cerr << "Error while opening file!" << endl;
}
