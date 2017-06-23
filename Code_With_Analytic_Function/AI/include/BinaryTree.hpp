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

class Node
{
    private:
        double m_key;
        bool m_isLeaf;
        int m_childType; // 0: fils gauche, 1: fils droite, -1: racine
        string m_code;

        Node* m_parent;
        Node* m_left;
        Node* m_right;

    public:
        Node(double key);
        Node(Node* node);
        ~Node() {};

        double key() { return m_key; };
        string code() { return m_code; };
        Node* parent() { return m_parent; };
        Node* left() { return m_left; };
        Node* right() { return m_right; };
        bool isLeaf() { return m_isLeaf; };
        int childType() { return m_childType; };

        void setKey(double key) { m_key = key; };
        void setCode(string code) { m_code = code; };
        void setParent(Node* node) { m_parent = node; };
        void setLeft(Node* node) { m_left = node; m_left->setChildType(0); };
        void setRight(Node* node) { m_right = node; m_right->setChildType(1); };
        void setIsLeaf(bool isLeaf) { m_isLeaf = isLeaf; };
        void setChildType(int t) { m_childType = t; };

        void displayNode();
        static void displayNodesRecursively(Node* node);
        static void clearNodesRecursively(Node* node);
        void clearNode();
};

class BinaryTree
{
    private:
        Node* m_root = NULL;
        vector<double> m_elems;

    public:
        BinaryTree();
        ~BinaryTree();

        Node* root() { return m_root; };
        vector<double>& elements() { return m_elems; };

        void addNode(double key);
        Node* searchNode(double key, double* key_inf, double* key_sup, bool lookAtChildren);
        double findKeySup(Node*, bool lookAtChildren);
        double findKeyInf(Node*, bool lookAtChildren);
        void clearTree();

        void displayBinaryTree();
        void tree2Vector(Node* node);

        static double getValueFromCode(string code);
        static string getParentCode(string code);
        static vector<string> computeChildrenCodes(string code);
};
typedef std::shared_ptr<BinaryTree> BinaryTreePtr;

#endif
