#ifndef BINARYTREE
#define BINARYTREE

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <memory>
#include <limits>
#include <iomanip>

using namespace std;

// Classe représentant un nœud
class Node
{
    private:
        double m_key; // Valeur du nœud
        bool m_isLeaf; // True si le nœud est une feuille
        int m_childType; // 0: fils gauche, 1: fils droite, -1: racine
        string m_code; // Code de Huffman du nœud

        Node* m_parent; // Nœud parent
        Node* m_left; // Fils gauche
        Node* m_right; // Fils droite

    public:
        Node(double key);
        Node(Node* node);
        ~Node() {};

	// Accesseurs
        double key() { return m_key; };
        string code() { return m_code; };
        Node* parent() { return m_parent; };
        Node* left() { return m_left; };
        Node* right() { return m_right; };
        bool isLeaf() { return m_isLeaf; };
        int childType() { return m_childType; };
	
	// mutateurs
        void setKey(double key) { m_key = key; };
        void setCode(string code) { m_code = code; };
        void setParent(Node* node) { m_parent = node; };
        void setLeft(Node* node) { m_left = node; m_left->setChildType(0); };
        void setRight(Node* node) { m_right = node; m_right->setChildType(1); };
        void setIsLeaf(bool isLeaf) { m_isLeaf = isLeaf; };
        void setChildType(int t) { m_childType = t; };

	// Fonctions d'affichage
        void displayNode();
        static void displayNodesRecursively(Node* node);

        static void clearNodesRecursively(Node* node);
        void clearNode();
};

// Classe implémentant un arbre binaire de recherche
// Utile pour la construction des points d'interpolation par dichotomie 
// Représente un ordre entre les points
// Un point p (nœud) est d'ordre plus grand qu'un point q (nœud) si q est un ancetre de p 
class BinaryTree
{
    private:
	// Racine de l'arbre (contient toujours la valeur 0.0) 
        Node* m_root = NULL;
	// Vecteur contenant les nœuds de l'arbre 
        vector<double> m_elems;

    public:
        BinaryTree();
        ~BinaryTree();

        Node* root() { return m_root; };
        vector<double>& elements() { return m_elems; };

	// Inserer l'élément de valeur key
        void addNode(double key);
	// Chercher le nœud de valeur key puis le nœud key_inf (resp. key_sup) qui correspond à la plus grande valeur inférieure 
	// (resp. plus petite valeur à key dont les valeurs 
        Node* searchNode(double key, double* key_inf, double* key_sup, bool lookAtChildren);
	// Chercher la plus petite valeur (nœud) supérieure à la valeur de node 
        double findKeySup(Node* node, bool lookAtChildren);
	// Chercher la plus grande valeur (nœud) inférieure à la valeur de node 
        double findKeyInf(Node* node, bool lookAtChildren);
        void clearTree();

	// Afficher l'arbre récursivement
        void displayBinaryTree();
        void tree2Vector(Node* node);
	
	// Trouver la valeur du nœud qui correspond au code Huffman code
        static double getValueFromCode(string code);
	// Trouver le code du nœud qui correspond au parent du nœud ayant code comme code de Huffman
        static string getParentCode(string code);
	// Calculer les codes des fils du nœud ayant code comme code de Huffman 
        static vector<string> computeChildrenCodes(string code);
	// Code de Huffman: soit le nœud n = DGDDG (chemin à partir de la racine) 
	// Le code de huffman correspondant est '10110' (1 si on va à droite, 0 si on va à gauche)

};
typedef std::shared_ptr<BinaryTree> BinaryTreePtr;

#endif
