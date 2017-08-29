#include "../include/BinaryTree.hpp"


/*********************************** Node *************************************/
Node::Node(double key)
{
    m_key = key;
    m_left = NULL;
    m_right = NULL;
    m_parent = NULL;
    m_isLeaf = true;
}
Node::Node(Node* node)
{
    m_key = node->key();
    m_left = node->left();
    m_right = node->right();
    m_parent = node->parent();
    m_isLeaf = node->isLeaf();
}
void Node::displayNode()
{
    // Affichage de toutes les informations dur un nœud
    cout << "[" << setprecision(numeric_limits<double>::digits10+1) << key();
    cout << ", ";
    if (parent()) cout << setprecision(numeric_limits<double>::digits10+1) << parent()->key();
    else cout << "X";
    cout << ", ";
    if (isLeaf()) cout << "X, X, True";
    else
    {
        if (left()) cout << setprecision(numeric_limits<double>::digits10+1) << left()->key();
        else cout << "X";
        cout << ", ";
        if (right()) cout << setprecision(numeric_limits<double>::digits10+1) << right()->key();
        else cout << "X";
        cout << ", ";
        cout << "False";
    }
    cout << ", ";
    if (m_childType == 1) cout << "D, ";
    else if (m_childType == 0) cout << "G, ";
    else cout << "Root" ;
    cout << code() << "]" << endl;
}
void Node::displayNodesRecursively(Node* node)
{
    // Affichage récursif d'un arbre
    if (node)
    {
        node->displayNode();
        displayNodesRecursively(node->left());
        displayNodesRecursively(node->right());
    }
}

void Node::clearNodesRecursively(Node* node)
{
    // Suppression récursive d'un arbre
    if (node!=NULL)
    {
        clearNodesRecursively(node->left());
        clearNodesRecursively(node->right());
        delete node;
    }
}

/******************************** BinaryTree **********************************/
BinaryTree::BinaryTree() {}
BinaryTree::~BinaryTree()
{
  clearTree();
}

void BinaryTree::addNode(double key)
{
  	// Insertion du nœud en fonction de sa valeur + calcul de son code de Huffman
  	// On commence à la racine
  	// si le nœud courant a une valeur plus grande on va à gauche
  	// sinon on va à droite
  	// Jusqu'à arriver à un noud null
  	Node *tmpNode;
    Node *tmpTree = m_root;

    // Nœud qu'on veut insérer
    Node *elem = new Node(key);
    // Initialisation du code du nouveau nœud
    string code = "";
    // si l'arbre n'est pas vide
    if (tmpTree)
    do
    {
        // Sauvegarde dans tmpNode l'adresse du dernier nœud visité
        tmpNode = tmpTree;
        // Mise à jour du parent du nouveau nœud
        elem->setParent(tmpNode);
        // Si la valeur du nouveau nœud est supérieure à la valeur du nœud courant
        if (key > tmpTree->key())
        {
            // on va à droite
            tmpTree = tmpTree->right();
            // Mise à jour du code du nouveau nœud
            code = code + "1";
            // Si feuille atteinte
            if (!tmpTree)
            {
                // Mise à jour du code est insertion à droite
                elem->setCode(code);
                tmpNode->setRight(elem);
                tmpNode->setIsLeaf(false);
            }
        }
        // Si la valeur du nouveau nœud est inférieure à la valeur du nœud courant
        else if (key < tmpTree->key())
        {
            // on va à gauche
            tmpTree = tmpTree->left();
            // Mise à jour du code du nouveau nœud
            code = code + "0";
            // Si feuille atteinte
            if (!tmpTree)
            {
                // Mise à jour du code est insertion à gauche
                elem->setCode(code);
                tmpNode->setLeft(elem);
                tmpNode->setIsLeaf(false);
            }
        }
        else
        {
            cout << "L'élément " << setprecision(numeric_limits<double>::digits10+1) << key << " existe déjà!" << endl;
            exit(1);
        }
    }
    while(tmpTree);
    else
    // si l'arbre est vide, on construit la racine
    {
        m_root = elem;
        m_root->setChildType(-1);
        m_root->setCode(code);
    }

}
double BinaryTree::getValueFromCode(string code)
{
    // Le code de 0.0 est la chaine vide ""
    if (code.compare("") == 0) return 0.0;
    // Le code de -1.0 est "0"
    if (code.compare("0") == 0) return -1.0;
    // Le code de 1.0 est "1"
    if (code.compare("1") == 0) return 1.0;
    double prev, cur, tmp;
    cur = (code[0]=='0') ? -1 : 1;
    prev = 0;
    // Voir la formule dans le rapport à la page 18
    for (int i=1; i<int(code.size()); i++)
    {
        tmp = cur;
        if (code[i]=='0') cur = cur - abs(prev-cur)/2;
        else cur = cur + abs(prev-cur)/2;
        prev = tmp;
    }
    return cur;
}
string BinaryTree::getParentCode(string code)
{
    code.pop_back();
    return code;
}

vector<string> BinaryTree::computeChildrenCodes(string code)
{
    // Code de Huffman: soit le nœud n = DGDDG (chemin à partir de la racine)
    // Le code de huffman correspondant est '10110' (1 si on va à droite, 0 si on va à gauche)
    // Le code de la racine est -1
    vector<string> childrenCodes;
    string codeLeft = code;
    string codeRight = code;
    codeLeft.push_back('0');
    codeRight.push_back('1');
    // Seul les nœuds de code "0" et "1" ont un childrenCodes de taille 1 car les valeurs descendantes de -1.0 (resp. 1.0)
    // sont uniquement -0.5 (resp. 0.5). Les autres nœuds ont toujours 2 fils possibles
    if (code.compare("0") == 0) childrenCodes.push_back(codeRight);
    else if (code.compare("1") == 0) childrenCodes.push_back(codeLeft);
    else
    {
        childrenCodes.push_back(codeLeft);
        childrenCodes.push_back(codeRight);
    }
    return childrenCodes;
}

Node* BinaryTree::searchNode(double key, double* key_inf, double* key_sup)
{
    Node *last_node = m_root, *temp = m_root;
    bool found = false;
	   // Boucle pour trouver le nœud key
    while(temp)
    {
        last_node = temp;
		    // si le nœud key est trouvé
        if (key == temp->key())
        {
			      // On cherche le nœuds inf et sup
			      // Pour avoir les valeurs des points les plus proches (qui entoure key)
            // But: Contruire la fonction de base (fonction chapeau par exemple)
            *key_sup = findKeySup(temp);
            *key_inf = findKeyInf(temp);
            found = true;
            return temp;
        }
    		// On va à droite ou à gauche selon la valeur du nœud courant
    		// A droite si inferieur à key
        if (key > temp->key())
            temp = temp->right();
    		// A gauche sinon
        else
            temp = temp->left();
    }
	  // Si le nœud n'est pas trouvé
    // On ne rentre jamais ici quand on veut construire une fonction de base (fonction chapeau par exemple)
    // car key correspond toujours à un point d'interpolation 1d présent naturellement dans l'arbre
    // preuve (dans la classe Interpolation et ses classe dérivées): appel depuis computeLastAlphaNu --> basisFunction_1D --> computeBoundariesForBasisFunction
    // le paramètre entré en premier lieu dans la fonction basisFunction_1D correspond à un point d'interpolation
    if (!found)
    {
		    // On cherche les valeurs voisines au dernier nœud visité (le plus proche à key)
        if (key<last_node->key())
        {
            *key_sup = last_node->key();
            *key_inf = findKeyInf(last_node);
        }
        else
        {
            *key_inf = last_node->key();
            *key_sup = findKeySup(last_node);
        }
    }
    return NULL;
}

double BinaryTree::findKeySup(Node* node)
{
    // On cherche le nœud ayant la plus petite valeur supérieure à node
    // On ne cherche pas parmi les descendants (voir computeBoundariesForBasisFunction() de la classe MixedInterpolation, PiecewiseInterpolation)

    // si node est la racine
    if (node->childType()==-1)
    {
        return node->key();
    }
    // s'il s'agit d'un fils gauche
    else if (node->childType()==0)
    {
        return node->parent()->key();
    }
    // s'il s'agit d'un fils droite
    else
    {
        Node* temp = node->parent();
        while (temp->childType()!=0)
        {
            if (temp->childType()==-1)
                return 1;
            temp = temp->parent();
        }
        return temp->parent()->key();
    }
}

double BinaryTree::findKeyInf(Node* node)
{
    // On cherche le nœud ayant la plus grande valeur inférieure à node
    // On ne cherche pas parmi les descendants (voir computeBoundariesForBasisFunction() de la classe MixedInterpolation, PiecewiseInterpolation)

    // si node est la racine
    if (node->childType()==-1)
    {
        return node->key();
    }
    // s'il s'agit d'un fils droite
    else if (node->childType()==1)
    {
        return node->parent()->key();
    }
    // s'il s'agit d'un fils gauche
    else
    {
        Node* temp = node->parent();
        while (temp->childType()!=1)
        {
            if (temp->childType()==-1)
                return -1;
            temp = temp->parent();
        }
        return temp->parent()->key();
    }
}

void BinaryTree::tree2Vector(Node* node)
{
    // Conversion de l'arbre en vecteur de réels
    if (node)
    {
        m_elems.push_back(node->key());
        tree2Vector(node->left());
        tree2Vector(node->right());
    }
}

void BinaryTree::clearTree()
{
  if (m_root)
  {
      Node::clearNodesRecursively(m_root);
      m_root = NULL;
  }
}

void BinaryTree::displayBinaryTree()
{
  if (m_root)
      Node::displayNodesRecursively(m_root);
}
