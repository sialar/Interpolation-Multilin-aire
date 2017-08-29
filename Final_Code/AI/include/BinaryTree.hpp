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

/**
 *  \file BinaryTree.hpp
 *  \brief Classe qui implémente un arbre binaire de recherche. Utile pour la version PiecewiseInterpolation
 *  \author SIALA Rafik
 *  \date 08/16
*/
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
        /**
          *  Construire un nœud à partir d'une valeur (une clé)
          *  \param key : valeur du nœud
        */
        Node(double key);
        /**
          *  Construire un nœud par copie
          *  \param node : nœud à copier
        */
        Node(Node* node);
        ~Node() {};

        /**
         *  Accesseur à l'attribut m_key
         *  \return valeur de m_key
        */
        double key() { return m_key; };
        /**
         *  Accesseur à l'attribut m_code
         *  \return valeur de m_code
        */
        string code() { return m_code; };
        /**
         *  Accesseur à l'attribut m_parent
         *  \return nœud parent
        */
        Node* parent() { return m_parent; };
        /**
         *  Accesseur à l'attribut m_left
         *  \return fils de gauche
        */
        Node* left() { return m_left; };
        /**
         *  Accesseur à l'attribut m_right
         *  \return fils de droite
        */
        Node* right() { return m_right; };
        /**
         *  Accesseur à l'attribut m_isLeff
         *  \return valeur de m_isLeaf
        */
        bool isLeaf() { return m_isLeaf; };
        /**
         *  Accesseur à l'attribut m_childType
         *  \return valeur de m_childType
        */
        int childType() { return m_childType; };

        /**
         *  Mutateur de l'attribut m_key
         *  \param key : valeur de la clé
        */
        void setKey(double key) { m_key = key; };
        /**
         *  Mutateur de l'attribut m_code
         *  \param code : code du nœud
        */
        void setCode(string code) { m_code = code; };
        /**
         *  Mutateur de l'attribut m_parent
         *  \param node : nœud parent
        */
        void setParent(Node* node) { m_parent = node; };
        /**
         *  Mutateur de l'attribut m_left
         *  \param node : fils de gauche
        */
        void setLeft(Node* node) { m_left = node; m_left->setChildType(0); };
        /**
         *  Mutateur de l'attribut m_right
         *  \param node : fils de droite
        */
        void setRight(Node* node) { m_right = node; m_right->setChildType(1); };
        /**
         *  Mutateur de l'attribut m_isLeaf
         *  \param isLeaf : valeur de m_isLeaf
        */
        void setIsLeaf(bool isLeaf) { m_isLeaf = isLeaf; };
        /**
         *  Mutateur de l'attribut m_childType
         *  \param t : valeur de m_childType
        */
        void setChildType(int t) { m_childType = t; };

        /**
         *  Afficher les informations sur le nœud
        */
        void displayNode();
        /**
         *  Afficher les informations sur tout les nœuds descendant de node
         *  d'une manière récursive
         *  \param node : racine d'un arbre
        */
        static void displayNodesRecursively(Node* node);

        /**
         *  Supprimer les nœuds descendant de node d'une manière récursive
         *  \param node : racine d'un arbre
        */
        static void clearNodesRecursively(Node* node);
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

        /**
         *  Accesseur à l'attribut m_root
         *  \return racine de l'arbre
        */
        Node* root() { return m_root; };
        /**
         *  Accesseur à l'attribut m_elems
         *  \return les elements de l'arbre
        */
        vector<double>& elements() { return m_elems; };

        /**
         *  Inserer l'élément de valeur key
         *  \param key : nouvelle valeur à inserer dans l'arbre
        */
        void addNode(double key);
        /**
         *  Chercher le nœud de valeur key puis le nœud key_inf (resp. key_sup) qui correspond à la plus grande valeur inférieure
         *  (resp. plus petite valeur supérieure) à key
         *  \param key : valeur du nœud cible
         *  \param key_inf : pointeur vers la plus grande valeur inférieure à key
         *  \param key_sup : pointeur vers la plus petite valeur supérieure à key
         *  \param return nœud de valeur key
        */
        Node* searchNode(double key, double* key_inf, double* key_sup);
        /**
         *  Chercher la plus petite valeur (nœud) supérieure à la valeur de node
         *  \param node : nœud concidéré
        */
        double findKeySup(Node* node);
        /**
         *  Chercher la plus grande valeur (nœud) inférieure à la valeur de node
         *  \param node : nœud concidéré
        */
        double findKeyInf(Node* node);
        /**
         *  Supprimer tous les éléments de l'arbre
        */
        void clearTree();

        /**
         *  Afficher l'arbre récursivement
        */
        void displayBinaryTree();
        /**
         *  Convertir l'arbre de racine node en vecteur de clés
         *  \param node : racine de l'arbre
        */
        void tree2Vector(Node* node);

        /**
         *  Trouver la valeur du nœud qui correspond au code de Huffman, code (voir page 18 du rapport)
         *  \param code : code de Huffman du nœud recherché
         *  \return : valeur du nœud recherché
        */
        static double getValueFromCode(string code);
        /**
         *  Trouver le code du nœud qui correspond au parent du nœud ayant code comme code de Huffman (voir page 18 du rapport)
         *  \param code : code de Huffman du fils du nœud recherché
         *  \return : code de Huffman du nœud parent
        */
        static string getParentCode(string code);
        /**
         *  Calculer les codes des fils du nœud ayant code comme code de Huffman (voir page 18 du rapport)
         *  \param code : code de Huffman du parent
         *  \return : vecteur (de taille <= 2) contenant les codes de Huffman des nœud fils
        */
        static vector<string> computeChildrenCodes(string code);
};

typedef std::shared_ptr<BinaryTree> BinaryTreePtr;

#endif
