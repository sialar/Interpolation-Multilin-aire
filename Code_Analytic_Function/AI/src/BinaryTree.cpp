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
    if (node)
    {
        node->displayNode();
        displayNodesRecursively(node->left());
        displayNodesRecursively(node->right());
    }
}

void Node::clearNodesRecursively(Node* node)
{
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
    Node *tmpNode;
    Node *tmpTree = m_root;

    Node *elem = new Node(key);
    string code = "";
    if (tmpTree)
    do
    {
        tmpNode = tmpTree;
        elem->setParent(tmpNode);
        if (key > tmpTree->key())
        {
            tmpTree = tmpTree->right();
            code = code + "1";
            if (!tmpTree)
            {
                elem->setCode(code);
                tmpNode->setRight(elem);
                tmpNode->setIsLeaf(false);
            }
        }
        else if (key < tmpTree->key())
        {
            tmpTree = tmpTree->left();
            code = code + "0";
            if (!tmpTree)
            {
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
    {
        m_root = elem;
        m_root->setChildType(-1);
        m_root->setCode(code);
    }

}
double BinaryTree::getValueFromCode(string code)
{
    if (code.compare("") == 0) return 0;
    if (code.compare("0") == 0) return -1;
    if (code.compare("1") == 0) return 1;
    double prev, cur, tmp;
    cur = (code[0]=='0') ? -1 : 1;
    prev = 0;
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
    vector<string> childrenCodes;
    string codeLeft = code;
    string codeRight = code;
    codeLeft.push_back('0');
    codeRight.push_back('1');
    if (code.compare("0") == 0) childrenCodes.push_back(codeRight);
    else if (code.compare("1") == 0) childrenCodes.push_back(codeLeft);
    else
    {
        childrenCodes.push_back(codeLeft);
        childrenCodes.push_back(codeRight);
    }
    return childrenCodes;
}

Node* BinaryTree::searchNode(double key, double* key_inf, double* key_sup, bool lookAtChildren)
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
