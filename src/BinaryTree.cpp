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
    Node* tmpNode = node;
    if (!node) return;
    if (tmpNode->left())  clearNodesRecursively(tmpNode->left());
    if (tmpNode->right()) clearNodesRecursively(tmpNode->right());
    delete tmpNode;
    node = NULL;
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
                m_hashmap.insert(pair<string,Node*>(code,elem));
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
                m_hashmap.insert(pair<string,Node*>(code,elem));
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
        m_root->setCode(code);
        m_hashmap.insert(pair<string,Node*>(code,m_root));
    }

}
double BinaryTree::getValueFromCode(string code)
{
    //map<string,Node*>::iterator it = m_hashmap.find(code);
    //if (it != m_hashmap.end()) return m_hashmap[code]->key();
    //else
    //{
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
    //}
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
    string code = "";
    while(temp)
    {
        last_node = temp;
        if (key == temp->key())
        {
            *key_sup = findKeySup(temp, lookAtChildren);
            *key_inf = findKeyInf(temp, lookAtChildren);
            //temp->setCode(code);
            found = true;
            return temp;
        }
        if (key > temp->key())
        {
            code = code + "1";
            temp = temp->right();
        }
        else
        {
            code = code + "0";
            temp = temp->left();
        }
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
  Node::clearNodesRecursively(m_root);
}

void BinaryTree::displayBinaryTree()
{
    Node::displayNodesRecursively(m_root);
}

int BinaryTree::getIndice(double l)
{
    if (l == 0) return 0;
    if (l == -1)  return 1;
    if (l == 1)  return 2;
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
double BinaryTree::getValue(int i)
{
  if (i == 0) return 0;
  if (i == 1)  return -1;
  if (i == 2)  return 1;

  int temp = i - 2; // nb points a partir de -0.5
  int e = 1;
  while (temp > pow(2,e))
  {
    temp -= pow(2,e);
    e++;
  }
  return (2*temp - 1 - pow(2,e)) / pow(2,e);
}
int BinaryTree::getParentIndice(int i)
{
    if (i < 3) return 0;
    if (i == 3) return 1;
    if (i == 4) return 2;

    int temp = i - 2;
    int e = 1;
    while (temp > pow(2,e))
    {
        temp -= pow(2,e);
        e++;
    }
    // le noeud i est le tempème noeud dans le niveau e
    // son parent est le numéro floor((temp+1)/2) dans le niveau (e-1)
    return (pow(2,e-1) + ((temp%2) ? (temp+1) : temp)/2);
}
double BinaryTree::getParentValue(int i)
{
    if (!i) return 2;
    return getValue(getParentIndice(i));
}
vector<double> BinaryTree::computeChildrenValue(int i)
{
    double parentVal, curentVal, offset;
    vector<double> childrenVal;
    parentVal = getParentValue(i);
    curentVal = getValue(i);
    offset = abs(parentVal - curentVal) / 2;
    if (i==1) childrenVal.push_back(curentVal + offset);
    else if (i==2) childrenVal.push_back(curentVal - offset);
    else
    {
        childrenVal.push_back(curentVal + offset);
        childrenVal.push_back(curentVal - offset);
    }
    return childrenVal;
}
bool BinaryTree::isAncestorOf(int i, int j) // i ancestor of j
{
    if (i>j) return false;
    if (i==j || i==0) return true;
    int temp = getParentIndice(j);
    while (temp!=0)
    {
        if (i==temp) return true;
        temp = getParentIndice(temp);
    }
    return false;
}
