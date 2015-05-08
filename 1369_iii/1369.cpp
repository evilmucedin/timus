#define _CRT_SECURE_NO_WARNINGS
#define M_PI 3.14159265358979323846

#include <cstdio>
#include <memory.h>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace std;

typedef vector<int> TIntVector;
typedef vector<double> TDoubles;

template<typename T>
T Sqr(T x)
{
    return x*x;
}

template<typename T>
T Max(T a, T b)
{
    return (a > b) ? a : b;
}

template<typename T>
T Min(T a, T b)
{
    return (a < b) ? a : b;
}

void Output(const TIntVector& vct)
{
    char buffer[100];
    for (size_t i = 0; i < vct.size(); ++i)
    {
        int num = vct[i] + 1;
        char* pBuffer = buffer;
        while (num)
        {
            *pBuffer = (num % 10) + '0';
            num /= 10;
            ++pBuffer;
        }
        *pBuffer = 0;
        std::reverse(buffer, pBuffer);
        *pBuffer = ' ';
        ++pBuffer;
        *pBuffer = 0;
        fputs(buffer, stdout);
    }
    fputs("\n", stdout);
}

typedef float TCoverNum;

class CoverTreePoint
{
private:
    TCoverNum _x;
    TCoverNum _y;
    int _index;

public:
    CoverTreePoint()
    {
    }

    CoverTreePoint(TCoverNum x, TCoverNum y, int index)
        : _x(x)
        , _y(y)
        , _index(index)
    {
    }

    CoverTreePoint(const CoverTreePoint& rhs)
        : _x(rhs._x)
        , _y(rhs._y)
        , _index(rhs._index)
    {
    }

    TCoverNum GetX() const
    {
        return _x;
    }

    TCoverNum GetY() const
    {
        return _y;
    }

    int GetIndex() const
    {
        return _index;
    }

    TCoverNum Distance(const CoverTreePoint& p) const;
    TCoverNum Distance2(const CoverTreePoint& p) const;
};

TCoverNum CoverTreePoint::Distance(const CoverTreePoint& p) const
{
    return sqrt(Sqr(_x - p._x) + Sqr(_y - p._y));
}

TCoverNum CoverTreePoint::Distance2(const CoverTreePoint& p) const
{
    return Sqr(_x - p._x) + Sqr(_y - p._y);
}

typedef vector<CoverTreePoint> CoverTreePoints;
    
class CoverTree
{
private:
    class CoverTreeNode
    {
    private:
        unordered_map< int, vector<CoverTreeNode*> > _childMap;
        CoverTreePoint _point;

    public:
        CoverTreeNode(const CoverTreePoint& p);
        void GetChildren(int level, std::vector<CoverTreeNode*>* result) const;
        const std::vector<CoverTreeNode*>& GetChildren(int level) const;
        void AddChild(int level, CoverTreeNode* p);
        
        const CoverTreePoint& GetPoint() const
        { 
            return _point;
        }
        
        TCoverNum Distance(const CoverTreeNode& p) const;

        bool HasPoint(const CoverTreePoint& p) const;

        vector<CoverTreeNode*> GetAllChildren() const;
    }; 

public:
    typedef vector<CoverTreeNode*> CoverTreeNodes;
    static const CoverTreeNodes EMPTY_NODES;

private:
    typedef std::pair<TCoverNum, CoverTreeNode*> TDistNodePair;

    CoverTreeNode* _root;
    unsigned int _numNodes;
    int _maxLevel;
    int _minLevel;

    CoverTreeNodes KNearestNodes(const CoverTreePoint& p, const unsigned int& k) const;
    CoverTreeNodes AllNearestNodes(const CoverTreePoint& p) const;
    bool InsertRec(const CoverTreePoint& p, const std::vector<TDistNodePair>& Qi, const int& level);

    TDistNodePair Distance(const CoverTreePoint& p, const std::vector<CoverTreeNode*>& Q);

public:
    static const TCoverNum base;

    CoverTree(const TCoverNum& maxDist, const CoverTreePoints& points); 
    ~CoverTree();

    void Insert(const CoverTreePoint& newPoint);

    CoverTreePoints KNearestNeighbors(const CoverTreePoint& p, const unsigned int& k) const;
    void AllNearestPoints(const CoverTreePoint& p, TIntVector* result) const;

    CoverTreeNode* GetRoot() const;
};

const TCoverNum CoverTree::base = 2.0;
const CoverTree::CoverTreeNodes CoverTree::EMPTY_NODES;

CoverTree::CoverTree(const TCoverNum& maxDist, const CoverTreePoints& points)
{
    _root = nullptr;
    _numNodes = 0;
    _maxLevel = static_cast<int>(ceil(log(maxDist)/log(base))) + 1;
    _minLevel = _maxLevel-1;
    for (CoverTreePoints::const_iterator it = points.begin(); it != points.end(); ++it)
    {
        Insert(*it);
    }
}

CoverTree::~CoverTree()
{
    if(!_root) return;
    std::vector<CoverTreeNode*> nodes;
    nodes.push_back(_root);
    while (!nodes.empty())
    {
        CoverTreeNode* byeNode = nodes[0];
        nodes.erase(nodes.begin());
        std::vector<CoverTreeNode*> children = byeNode->GetAllChildren();
        nodes.insert(nodes.begin(), children.begin(), children.end());
        delete byeNode;
    }   
}

CoverTree::CoverTreeNodes CoverTree::KNearestNodes(const CoverTreePoint& p, const unsigned int& k) const
{
    if (!_root)
    {
        return CoverTreeNodes();
    }
    TCoverNum maxDist = p.Distance(_root->GetPoint());
    set<TDistNodePair> minNodes;

    minNodes.insert(make_pair(maxDist,_root));
    std::vector<TDistNodePair> Qj(1, make_pair(maxDist,_root));
    for(int level = _maxLevel; level>=_minLevel; level--) {
        std::vector<TDistNodePair>::const_iterator it;
        int size = Qj.size();
        for(int i=0; i<size; i++) {
            std::vector<CoverTreeNode*> children = Qj[i].second->GetChildren(level);
            for(CoverTreeNodes::const_iterator it2=children.begin(); it2!=children.end(); ++it2) {
                TCoverNum d = p.Distance((*it2)->GetPoint());
                if(d < maxDist || minNodes.size() < k)
                {
                    minNodes.insert(make_pair(d,*it2));
                    if (minNodes.size() > k)
                    {
                        minNodes.erase(--minNodes.end());
                    }
                    maxDist = (--minNodes.end())->first;
                }
                Qj.push_back(make_pair(d,*it2));
            }
        }
        TCoverNum sep = maxDist + pow(base, level);
        size = Qj.size();
        for(int i=0; i<size; i++) {
            if(Qj[i].first > sep) {
                //quickly removes an element from a vector w/o preserving order.
                Qj[i]=Qj.back();
                Qj.pop_back();
                size--; i--;
            }
        }
    }
    
    CoverTreeNodes kNN;
    for (set<TDistNodePair>::const_iterator it = minNodes.begin(); it != minNodes.end(); ++it)
    {
        kNN.push_back(it->second);
    }
    return kNN;
}

static const TCoverNum EPS = TCoverNum(1e-8);

CoverTree::CoverTreeNodes CoverTree::AllNearestNodes(const CoverTreePoint& p) const
{
    if (!_root)
    {
        return CoverTreeNodes();
    }
    TCoverNum maxDist = p.Distance2(_root->GetPoint());
    CoverTreeNodes minNodes;
    minNodes.push_back(_root);
    vector<TDistNodePair> q(1, make_pair(maxDist,_root));
    for(int level = _maxLevel; level >= _minLevel; --level)
    {
        int size = q.size();
        for (int i = 0; i < size; ++i)
        {
            const CoverTreeNodes& children = q[i].second->GetChildren(level);
            for (CoverTreeNodes::const_iterator it2 = children.begin(); it2 != children.end(); ++it2)
            {
                TCoverNum d = p.Distance2((*it2)->GetPoint());
                if (d < maxDist + EPS)
                {
                    minNodes.push_back(*it2);
                    if (d < maxDist)
                    {
                        maxDist = d;
                    }
                }
                q.push_back(make_pair(d, *it2));
            }
        }
        const TCoverNum sep = Sqr(sqrt(maxDist) + pow(base, level)) + EPS;
        size = q.size();
        for (int i = 0; i < size; ++i)
        {
            if (q[i].first > sep)
            {
                q[i] = q.back();
                q.pop_back();
                --size;
                --i;
            }
        }
    }
    
    return minNodes;
}

static const TCoverNum INF = 1e6f;

bool CoverTree::InsertRec(const CoverTreePoint& p, const std::vector<TDistNodePair>& Qi, const int& level)
{
    std::vector<std::pair<TCoverNum, CoverTreeNode*> > Qj;
    TCoverNum sep = pow(base, level);
    double minDist = INF;
    std::pair<double, CoverTreeNode*> minQiDist(minDist, nullptr);
    for(std::vector< std::pair<TCoverNum, CoverTreeNode*> >::const_iterator it = Qi.begin(); it != Qi.end(); ++it)
    {
        if(it->first < minQiDist.first)
            minQiDist = *it;
        if(it->first<minDist)
            minDist = it->first;
        if(it->first <= sep)
            Qj.push_back(*it);
        CoverTreeNodes children = it->second->GetChildren(level);
        for(CoverTreeNodes::const_iterator it2 = children.begin(); it2 != children.end(); ++it2)
        {
            double d = p.Distance((*it2)->GetPoint());
            if (d < minDist)
            {
                minDist = d;
            }
            if (d <= sep)
            {
                Qj.push_back(make_pair(d,*it2));
            }
        }
    }
    if (minDist > sep)
    {
        return true;
    }
    else
    {
        bool found = InsertRec(p, Qj, level-1);
        //distNodePair minQiDist = distance(p,Qi);
        if (found && minQiDist.first <= sep)
        {
            if (level < _minLevel + 1) 
            {
                _minLevel = level-1;
            }
            minQiDist.second->AddChild(level, new CoverTreeNode(p));
            ++_numNodes;
            return false;
        }
        else 
        {
            return found;
        }
    }
}

std::pair<TCoverNum, CoverTree::CoverTreeNode*> CoverTree::Distance(const CoverTreePoint& p, const std::vector<CoverTreeNode*>& Q)
{
    double minDist = INF;
    CoverTreeNode* minNode;
    for (CoverTreeNodes::const_iterator it=Q.begin();it!=Q.end(); ++it)
    {
        double dist = p.Distance((*it)->GetPoint());
        if (dist < minDist)
        {
            minDist = dist;
            minNode = *it;
        }
    }
    return make_pair(minDist,minNode);  
}

void CoverTree::Insert(const CoverTreePoint& newPoint)
{
    if (!_root)
    {
        _root = new CoverTreeNode(newPoint);
        _numNodes = 1;
    }
    else
    {
        InsertRec(newPoint, std::vector<TDistNodePair>(1, make_pair(_root->Distance(newPoint), _root)), _maxLevel);           
    }
}

CoverTreePoints CoverTree::KNearestNeighbors(const CoverTreePoint& p, const unsigned int& k) const
{
    if (_root == NULL)
        return CoverTreePoints();
    CoverTreeNodes v = KNearestNodes(p, k);
    CoverTreePoints kNN;
    for (CoverTreeNodes::const_iterator it=v.begin(); it!=v.end(); ++it)
    {
        const CoverTreePoint& p = (*it)->GetPoint();
        kNN.push_back(p);
        if (kNN.size() >= k)
        {
            break;
        }
    }
    return kNN;
}

void CoverTree::AllNearestPoints(const CoverTreePoint& p, TIntVector* result) const
{
    CoverTreeNodes v = AllNearestNodes(p);
    result->resize(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        (*result)[i] = v[i]->GetPoint().GetIndex();
    }
}

CoverTree::CoverTreeNode* CoverTree::GetRoot() const
{
    return _root;
}

CoverTree::CoverTreeNode::CoverTreeNode(const CoverTreePoint& p)
    : _point(p)
{
}

const CoverTree::CoverTreeNodes& CoverTree::CoverTreeNode::GetChildren(int level) const
{
    unordered_map<int, vector<CoverTreeNode*> >::const_iterator it = _childMap.find(level);
    if (it != _childMap.end())
    {
        return it->second;
    }
    return EMPTY_NODES;
}

void CoverTree::CoverTreeNode::AddChild(int level, CoverTreeNode* p)
{
    _childMap[level].push_back(p);
}

TCoverNum CoverTree::CoverTreeNode::Distance(const CoverTreeNode& p) const
{
    return _point.Distance(p.GetPoint());
}

CoverTree::CoverTreeNodes CoverTree::CoverTreeNode::GetAllChildren() const
{
    CoverTreeNodes children;
    for (unordered_map<int, CoverTreeNodes>::const_iterator it = _childMap.begin(); it!=_childMap.end(); ++it)
    {
        children.insert(children.end(), it->second.begin(), it->second.end());
    }
    return children;
}


void GenBig()
{
    FILE* fOut = fopen("big.txt", "w");
    static const int N = 100000;
    fprintf(fOut, "%d\n", N);
    for (int i = 0; i < N; ++i)
    {
        long double ldi = i;
        long double angle = ldi/N*2.0*M_PI;
        static const long double R = 10000.0;
        long double x = R*cos(angle);
        long double y = R*sin(angle);
        fprintf(fOut, "%Lf %Lf\n", x, y);
    }
    static const int M = 10000;
    fprintf(fOut, "%d\n", M);
    for (int i = 0; i < M; ++i)
    {
        long double ldi = i;
        long double angle = ldi/M*2.0*M_PI;
        static const long double R = 0.1;
        long double x = R*cos(angle);
        long double y = R*sin(angle);
        fprintf(fOut, "%Lf %Lf\n", x, y);
    }
    fclose(fOut);
}

void GenInt()
{
    FILE* fOut = fopen("int.txt", "w");
    static const int N = 10000;
    fprintf(fOut, "%d\n", N);
    for (int i = 0; i < N; ++i)
    {
        fprintf(fOut, "%d %d\n", rand() % 10, rand() % 10);
    }
    static const int M = 10000;
    fprintf(fOut, "%d\n", M);
    for (int i = 0; i < M; ++i)
    {
        fprintf(fOut, "%Lf %Lf\n", rand() % 10, rand() % 10);
    }
    fclose(fOut);
}

int main()
{
#ifndef ONLINE_JUDGE
    // GenBig();
    // freopen("big.txt", "r", stdin);
    freopen("input.txt", "r", stdin);
    // freopen("small.txt", "r", stdin);
    // GenInt();
    // freopen("int.txt", "r", stdin);
#endif

    int m;
    scanf("%d", &m);

    int size = m + 8;

    long double* xcd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    long double* ycd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    long double mx = 1.0;
    for (int i = 0; i < m; ++i)
    {
        scanf("%Lf%Lf", &xcd[i], &ycd[i]);
        mx = Max(mx, abs(xcd[i]));
        mx = Max(mx, abs(ycd[i]));
    }

    TIntVector indices(m);
    for (int i = 0; i < m; ++i)
    {
        indices[i] = i;
    }

    for (int i = 0; i < m; ++i)
    {
        swap(indices[i], indices[i + (rand() % (m - i))]);
    }

    TDoubles xc(m);
    TDoubles yc(m);
    for (int i = 0; i < m; ++i)
    {
        xc[i] = xcd[indices[i]] / mx;
        yc[i] = ycd[indices[i]] / mx;
    }

    CoverTreePoints points(m);
    for (int i = 0; i < m; ++i)
    {
        points[i] = CoverTreePoint(xc[i], yc[i], indices[i]);
    }

    CoverTree tree(10, points);

    int n;
    scanf("%d", &n);
    vector<int> result;
    result.reserve(m);
    vector<int> result2;
    result2.reserve(m);
    for (int i = 0; i < n; ++i)
    {
        long double xd, yd;
        scanf("%Lf%Lf", &xd, &yd);

        double x = xd/mx;
        double y = yd/mx;

        CoverTreePoint q(x, y, -1);

        tree.AllNearestPoints(q, &result);

        result2.clear();
        long double mind = 1e15;
        long double mindMin = mind;
        long double mindMax = mind;
        for (int j = 0; j < result.size(); ++j)
        {
            int realindex = result[j];
            long double dist = Sqr(xd - xcd[realindex]) + Sqr(yd - ycd[realindex]);
            static const long double LDEPS = 1e-10;
            if (dist < mindMin)
            {
                mind = dist;
                mindMin = dist - LDEPS;
                mindMax = dist + LDEPS;
                result2.clear();
            }
            if (dist <= mindMax)
            {
                result2.push_back(realindex);
            }
        }

        std::sort(result2.begin(), result2.end());
        result2.erase(unique(result2.begin(), result2.end()), result2.end());

        Output(result2);
    }

    return 0;
}