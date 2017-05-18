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
    cout << "[" << key();
    cout << ", ";
    if (parent()) cout << parent()->key();
    else cout << "X";
    cout << ", ";
    if (isLeaf()) cout << "X, X, True";
    else
    {
        if (left()) cout << left()->key();
        else cout << "X";
        cout << ", ";
        if (right()) cout << right()->key();
        else cout << "X";
        cout << ", ";
        cout << "False";
    }
    cout << ", ";
    if (m_childType == 1) cout << "D]" << endl;
    else if (m_childType == 0) cout << "G]" << endl;
    else cout << "Root]" << endl;
}
void Node::displayNodesRecursively(Node* node)
{
    if (node)
    {
        node->displayNode();
        displayNodesRecursively(node->left());
        displayNodesRecursively(node->right());
    }
}


/******************************** BinaryTree **********************************/
BinaryTree::BinaryTree()
{
    addNode(0);
}

void BinaryTree::initTree(int depth)
{
    double key, denom;
    for (int i=0; i<=depth; i++)
    {
        denom = pow(2,i);
        for (int j=0; j<denom; j++)
        {
            if (j%2)
            {
                key = j/denom;
                addNode(key);
                addNode(-key);
            }
        }
    }
}

void BinaryTree::addNode(double key)
{
    Node *tmpNode;
    Node *tmpTree = m_root;

    Node *elem = new Node(key);

    if (tmpTree)
    do
    {
        tmpNode = tmpTree;
        elem->setParent(tmpNode);
        if (key > tmpTree->key())
        {
            tmpTree = tmpTree->right();
            if (!tmpTree)
            {
                tmpNode->setRight(elem);
                tmpNode->setIsLeaf(false);
            }
        }
        else if (key < tmpTree->key())
        {
            tmpTree = tmpTree->left();
            if (!tmpTree)
            {
                tmpNode->setLeft(elem);
                tmpNode->setIsLeaf(false);
            }
        }
        else
            cout << "L'élément " << key << " existe déjà!" << endl;
    }
    while(tmpTree);
    else
    {
        m_root = elem;
        m_root->setChildType(-1);
    }

}

Node* BinaryTree::searchNode(double key, double* key_sup, double* key_inf, bool lookAtChildren)
{
    Node *last_node = m_root, *temp = m_root;
    bool found = false;
    while(temp)
    {
        last_node = temp;
        if (key == temp->key())
        {
            *key_sup = findKeySup(temp, lookAtChildren);
            *key_inf = findKeyInf(temp, lookAtChildren);
            found = true;
            return temp;
        }
        if (key > temp->key())
            temp = temp->right();
        else
            temp = temp->left();
    }
    if (!found)
    {
        if (key<last_node->key())
        {
            *key_sup = last_node->key();
            *key_inf = findKeyInf(last_node, lookAtChildren);
        }
        else
        {
            *key_inf = last_node->key();
            *key_sup = findKeySup(last_node, lookAtChildren);
        }
    }
    return NULL;
}


double BinaryTree::findKeySup(Node* node, bool lookAtChildren)
{
    if (node->right() && lookAtChildren)
    {
        Node* temp = node->right();
        while (temp->left())
            temp = temp->left();
        return temp->key();
    }
    if (node->childType()==-1)
    {
        return node->key();
    }
    else if (node->childType()==0)
    {
        return node->parent()->key();
    }
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

double BinaryTree::findKeyInf(Node* node, bool lookAtChildren)
{
    if (node->left() && lookAtChildren)
    {
        Node* temp = node->left();
        while (temp->right())
            temp = temp->right();
        return temp->key();
    }
    if (node->childType()==-1)
    {
        return node->key();
    }
    else if (node->childType()==1)
    {
        return node->parent()->key();
    }
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
    if (node)
    {
        m_elems.push_back(node->key());
        tree2Vector(node->left());
        tree2Vector(node->right());
    }
}

void BinaryTree::displayBinaryTree()
{
    Node::displayNodesRecursively(m_root);
}

int BinaryTree::getIndice(double l)
{
    if (l == -1) return 0;
    if (l == 1)  return 1;
    if (l == 0)  return 2;
    double temp = l;
    int n = 0;
    while (floor(temp) != temp)
    {
        temp *= 2;
        n++;
    }
    int start = 2;
    for (int k=1; k<n; k++)
        start += pow(2,n-k);
    double a = 0.5 * (temp + pow(2,n)) + 1;
    return start + int(a);
}
