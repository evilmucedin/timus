#define _CRT_SECURE_NO_WARNINGS

#include <immintrin.h>
#include <intrin.h>

#include <cstdio>
#include <cmath>
#include <cassert>

#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <memory>
#include <iostream>
#include <list>

#define M_PI 3.14159265358979323846

using namespace std;

static const double INF = 1e9;

typedef vector<int> Integers;
typedef set<int> IntegerSet;

class Exception : public std::exception
{
private:
    const char* _s;

public:
    Exception(const char* s)
        : _s(s)
    {
    }

    const char* what() const throw() {
        return _s;
    }
};

template<typename T>
T Sqr(T x) {
    return x*x;
}

void Output(const Integers& vct) {
    char buffer[100];
    for (size_t i = 0; i < vct.size(); ++i) {
        int num = vct[i] + 1;
        char* pBuffer = buffer;
        while (num) {
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

typedef double TBase;

struct Vector;

struct Point
{
    TBase _x;
    TBase _y;
    int _index;

    Point() {

    }

    Point(TBase x, TBase y, int index)
        : _x(x)
        , _y(y)
        , _index(index)
    {

    }

    TBase Distance2(const Point& p) {
        return Sqr(_x - p._x) + Sqr(_y - p._y);
    }
};

typedef vector<Point> Points;

struct Vector
{
    TBase _dx;
    TBase _dy;

    Vector(TBase dx, TBase dy) 
        : _dx(dx)
        , _dy(dy)
    {

    }

    TBase Length2() const {
        return Sqr(_dx) + Sqr(_dy);
    }

    TBase Length() const {
        return sqrt(Length2());
    }

    void Normalize() {
        double norm = Length();
        _dx /= norm;
        _dy /= norm;
    }
};

TBase Dot(const Point& p, const Vector& v) {
    return p._x*v._dx + p._y*v._dy;
}

TBase Angle(const Point& p, const Point& p0) {
    return atan2(p._y - p0._y, p._x - p0._x);
}

struct PointProjection
{
    const Point* _p;
    TBase _projection;

    PointProjection() {
    }

    PointProjection(const Point* p, TBase projection)
        : _p(p)
        , _projection(projection)
    {

    }

    bool operator<(const PointProjection& rhs) const
    {
        return _projection < rhs._projection;
    }
};
typedef vector<PointProjection> PointProjections;

void GenBig()
{
    FILE* fOut = fopen("big.txt", "w");
    static const int N = 100000;
    fprintf(fOut, "%d\n", N);
    for (int i = 0; i < N; ++i)
    {
        long double ldi = i;
        long double angle = ldi / N*2.0*M_PI;
        static const long double R = 10000.0;
        long double x = R*cos(angle);
        long double y = R*sin(angle);
        fprintf(fOut, "%Lf %Lf\n", x, y);
    }
    static const int M = 10000;
    fprintf(fOut, "%d\n", M);
    for (int i = 0; i < M; ++i) {
        long double ldi = i;
        long double angle = ldi / M*2.0*M_PI;
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

struct DPoint {
    long double _x;
    long double _y;

    DPoint()
    {

    }

    DPoint(long double x, long double y)
        : _x(x)
        , _y(y)
    {

    }

    long double distance2(const DPoint& p) const {
        return Sqr(p._x - _x) + Sqr(p._y - _y);
    }
};

typedef vector<DPoint> DPoints;

double Rand() {
    return 2.0*(static_cast<double>(rand()) / RAND_MAX - 0.5);
}

/*
void DedupGraph(vector<Integers>& graph) {
    for (int i = 0; i < graph.size(); ++i) {
        sort(graph[i].begin(), graph[i].end());
        graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());
    }
}
*/

int main() {
#ifndef ONLINE_JUDGE
    GenBig();
    freopen("big.txt", "r", stdin);
    // freopen("input.txt", "r", stdin);
#endif

    // cout << sizeof(Triangle) << endl;

    int m;
    scanf("%d", &m);
    DPoints dpoints(m);
    Integers indices(m);
    for (int i = 0; i < m; ++i) {
        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        dpoints[i] = DPoint(x, y);
        indices[i] = i;
    }

    for (int i = 0; i < m; ++i) {
        swap(indices[i], indices[i + (rand() % (m - i))]);
    }

    vector<IntegerSet> graph;
    graph.resize(m);
    vector<bool> used;
    used.resize(m);
    queue<int> qu;

    Points points(m);
    for (int i = 0; i < m; ++i) {
        points[i] = Point(dpoints[i]._x, dpoints[i]._y, i);
    }
   
    PointProjections projections(m);
    for (int i = 0; i < m; ++i) {
        projections[i]._p = &points[i];
    }

    static const int D = 1;

    for (int it = 0; it < 300; ++it) {
        Vector vRand(Rand(), Rand());
        vRand.Normalize();
        for (int i = 0; i < m; ++i) {
            projections[i]._projection = Dot(points[i], vRand);
        }
        sort(projections.begin(), projections.end());
        for (int i = 0; i < m; ++i) {
            for (int di = -D; di <= D; ++di) {
                if (di) {
                    int nextIndex = i + di;
                    while (nextIndex < m) {
                        nextIndex += m;
                    }
                    while (nextIndex >= m) {
                        nextIndex -= m;
                    }
                    graph[projections[i]._p->_index].insert(projections[nextIndex]._p->_index);
                }
            }
        }
    }

    // DedupGraph(graph);

    for (int it = 0; it < 300; ++it) {
        int index = rand() % m;
        for (int i = 0; i < m; ++i) {
            projections[i]._projection = Angle(points[i], points[index]);
        }
        sort(projections.begin(), projections.end());
        for (int i = 0; i < m; ++i) {
            for (int di = -D; di <= D; ++di) {
                if (di) {
                    int nextIndex = i + di;
                    while (nextIndex < m) {
                        nextIndex += m;
                    }
                    while (nextIndex >= m) {
                        nextIndex -= m;
                    }
                    graph[projections[i]._p->_index].insert(projections[nextIndex]._p->_index);
                }
            }
        }
    }

    // DedupGraph(graph);

    int n;
    scanf("%d", &n);
    Integers result;
    Integers dindices;
    Integers temp;

    for (int i = 0; i < n; ++i) {        
        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        DPoint dq(x, y);

        dindices.clear();
        temp.clear();
        TBase mind = INF;
        for (int i = 0; i < 10; ++i) {
            int index = rand() % m;
            if (!used[index]) {
                qu.push(index);
                used[index] = true;
                temp.push_back(index);
            }
        }
        while (!qu.empty()) {
            int index = qu.front();
            qu.pop();
            TBase dist = dpoints[index].distance2(dq);
            if (dist <= mind + 1e-9) {
                if (dist < mind) {
                    if (dist + 1e-9 < mind) {
                        result.clear();
                    }
                    mind = dist;
                }
                dindices.push_back(index);
                for (IntegerSet::const_iterator toItem = graph[index].begin(); toItem != graph[index].end(); ++toItem) {
                    int v = *toItem;
                    if (!used[v]) {
                        temp.push_back(v);
                        used[v] = true;
                        qu.push(v);
                    }
                }
            }
        }
        for (int k = 0; k < temp.size(); ++k) {
            used[temp[k]] = false;
        }
        sort(dindices.begin(), dindices.end());

        // BigDecimal min = new BigDecimal(1e12);
        result.clear();
        long double min = 1e12;
        int prevIndex = -1;
        for (int k = 0; k < dindices.size(); ++k) {
            int index = dindices[k];
            if (index > prevIndex) {
                long double d = dpoints[index].distance2(dq);
                if (d < min) {
                    min = d;
                    result.clear();
                }
                if (d == min) {
                    result.push_back(index);
                }
                prevIndex = index;
            }
        }
 
        Output(result);
    }

    return 0;
}