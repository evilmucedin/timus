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

union TVector
{
    __m128 v;
    float f[4];
};

typedef vector<float> TFloats;

static const float FINF = 1e8f;

static void ConvertVector(const TFloats& floats, TVector** result, int* len)
{
    const int m = static_cast<int>(floats.size());
    *len = m / 4;
    if (m % 4)
    {
        ++(*len);
    }
    *result = (TVector*)_mm_malloc(sizeof(TVector)*(*len), 32);
    for (int i = 0; i < *len; ++i)
    {
        float values[] = { FINF, FINF, FINF, FINF };
        for (int j = 0; j < 4; ++j)
        {
            if (4 * i + j < m)
            {
                values[j] = floats[4 * i + j];
            }
        }
        (*result)[i].v = _mm_load_ps(values);
    }
}

static double Median(const TFloats& floats)
{
    TFloats temp(floats);
    sort(temp.begin(), temp.end());
    return temp[temp.size() / 2];
}

struct KDTree
{
    int _len;
    TVector* _x;
    TVector* _y;
    bool _isLeaf;
    KDTree* _left;
    KDTree* _right;
    Integers _indices;
    float _minX;
    float _maxX;
    float _minY;
    float _maxY;

    KDTree(int depth, const TFloats& x, const TFloats& y, const Integers& indices)
        : _left(nullptr)
        , _right(nullptr)
    {
        if (depth == 9)
        {
            ConvertVector(x, &_x, &_len);
            ConvertVector(y, &_y, &_len);
            _indices = indices;
            _isLeaf = true;
            _minX = FINF;
            _maxX = -FINF;
            _minY = FINF;
            _maxY = -FINF;
            for (size_t i = 0; i < x.size(); ++i)
            {
                _minX = Min(_minX, x[i]);
                _maxX = Max(_maxX, x[i]);
                _minY = Min(_minY, y[i]);
                _maxY = Max(_maxY, y[i]);
            }
        }
        else
        {
            float separator = (depth & 1) ? Median(x) : Median(y);

            TFloats xLeft;
            TFloats yLeft;
            Integers indicesLeft;
            TFloats xRight;
            TFloats yRight;
            Integers indicesRight;

            for (size_t i = 0; i < x.size(); ++i)
            {
                bool toLeft = (depth & 1) ? (x[i] < separator) : (y[i] < separator);
                if (toLeft)
                {
                    xLeft.push_back(x[i]);
                    yLeft.push_back(y[i]);
                    indicesLeft.push_back(indices[i]);
                }
                else
                {
                    xRight.push_back(x[i]);
                    yRight.push_back(y[i]);
                    indicesRight.push_back(indices[i]);
                }
            }

            if (xLeft.size())
            {
                _left = new KDTree(depth + 1, xLeft, yLeft, indicesLeft);
            }
            if (xRight.size())
            {
                _right = new KDTree(depth + 1, xRight, yRight, indicesRight);
            }

            _isLeaf = false;
        }
    }

#ifndef _MSC_VER
    __attribute__((force_align_arg_pointer))
#endif
        void Solve(float x, float y, const TVector& x4, const TVector& y4, float* minDist, TVector* minMax, float* minMin, Integers* result, int* limit) const
    {
            if (*limit < 0)
            {
                return;
            }

            if (_isLeaf)
            {
                if (!_len)
                {
                    return;
                }

                float d;
                if (x < _minX)
                {
                    if (y < _minY)
                    {
                        d = Min(_minX - x, _minY - y);
                    }
                    else if (y > _maxY)
                    {
                        d = Min(_minX - x, y - _maxY);
                    }
                    else
                    {
                        d = _minX - x;
                    }
                }
                else if (x > _maxX)
                {
                    if (y < _minY)
                    {
                        d = Min(x - _maxX, _minY - y);
                    }
                    else if (y > _maxY)
                    {
                        d = Min(x - _maxX, y - _maxY);
                    }
                    else
                    {
                        d = x - _maxX;
                    }
                }
                else
                {
                    if (y < _minY)
                    {
                        d = _minY - y;
                    }
                    else if (y > _maxY)
                    {
                        d = y - _maxY;
                    }
                    else
                    {
                        d = 0;
                    }
                }

                if (d*d > *minDist)
                {
                    return;
                }

                static const float EPS = 4e-7f;
                *limit -= _len;
                for (int j = 0; j < _len; ++j)
                {
                    TVector dx4;
                    dx4.v = _mm_sub_ps(_x[j].v, x4.v);
                    dx4.v = _mm_mul_ps(dx4.v, dx4.v);
                    TVector dy4;
                    dy4.v = _mm_sub_ps(_y[j].v, y4.v);
                    dy4.v = _mm_mul_ps(dy4.v, dy4.v);
                    dx4.v = _mm_add_ps(dx4.v, dy4.v);

                    TVector cmpMax;
                    cmpMax.v = _mm_cmplt_ps(dx4.v, minMax->v);
                    for (int k = 0; k < 4; ++k)
                    {
                        if (cmpMax.f[k])
                        {
                            if (dx4.f[k] < *minDist)
                            {
                                if (dx4.f[k] < *minMin)
                                {
                                    result->clear();
                                }
                                *minDist = dx4.f[k];
                                minMax->v = _mm_set1_ps(*minDist + EPS);
                                *minMin = *minDist - EPS;
                            }
                            result->push_back(_indices[4 * j + k]);
                        }
                    }
                }
            }
            else
            {
                if (_left)
                {
                    _left->Solve(x, y, x4, y4, minDist, minMax, minMin, result, limit);
                }
                if (_right)
                {
                    _right->Solve(x, y, x4, y4, minDist, minMax, minMin, result, limit);
                }
            }
        }
};

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

void DedupGraph(vector<Integers>& graph) {
    for (int i = 0; i < graph.size(); ++i) {
        sort(graph[i].begin(), graph[i].end());
        graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());
    }
}

int main() {
#ifndef ONLINE_JUDGE
    // GenBig();
    // freopen("big.txt", "r", stdin);
    freopen("input.txt", "r", stdin);
#endif

    // cout << sizeof(Triangle) << endl;

    int m;
    scanf("%d", &m);
    DPoints dpoints(m);
    long double mx = 1e-6;
    for (int i = 0; i < m; ++i) {
        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        dpoints[i] = DPoint(x, y);
        mx = Max(mx, abs(x));
        mx = Max(mx, abs(y));
    }

    vector<Integers> graph;
    graph.resize(m);

    {
        Points points(m);
        for (int i = 0; i < m; ++i) {
            points[i] = Point(dpoints[i]._x, dpoints[i]._y, i);
        }

        PointProjections projections(m);
        for (int i = 0; i < m; ++i) {
            projections[i]._p = &points[i];
        }

        static const int N_IT = 7;
        static const int D = 1;

        for (int it = 0; it < N_IT; ++it) {
            Vector vRand(Rand(), Rand());
            vRand.Normalize();
            for (int i = 0; i < m; ++i) {
                projections[i]._projection = Dot(points[i], vRand);
            }
            sort(projections.begin(), projections.end());
            for (int i = 0; i < m; ++i) {
                for (int di = 1; di <= D; ++di) {
                    int nextIndex = i + di;
                    while (nextIndex < m) {
                        nextIndex += m;
                    }
                    while (nextIndex >= m) {
                        nextIndex -= m;
                    }
                    int index1 = projections[i]._p->_index;
                    int index2 = projections[nextIndex]._p->_index;
                    graph[index1].push_back(index2);
                    graph[index2].push_back(index1);
                }
            }
        }

        DedupGraph(graph);

        for (int it = 0; it < N_IT; ++it) {
            int index = rand() % m;
            for (int i = 0; i < m; ++i) {
                projections[i]._projection = Angle(points[i], points[index]);
            }
            sort(projections.begin(), projections.end());
            for (int i = 0; i < m; ++i) {
                for (int di = 1; di <= D; ++di) {
                    int nextIndex = i + di;
                    while (nextIndex < m) {
                        nextIndex += m;
                    }
                    while (nextIndex >= m) {
                        nextIndex -= m;
                    }
                    int index1 = projections[i]._p->_index;
                    int index2 = projections[nextIndex]._p->_index;
                    graph[index1].push_back(index2);
                    graph[index2].push_back(index1);
                }
            }
        }

        DedupGraph(graph);
    }

    KDTree* kdTree;
    {
        Integers indices(m);
        for (int i = 0; i < m; ++i) {
            indices[i] = i;
        }
        for (int i = 0; i < m; ++i) {
            swap(indices[i], indices[i + (rand() % (m - i))]);
        }

        TFloats xc(m);
        TFloats yc(m);
        for (int i = 0; i < m; ++i)
        {
            xc[i] = dpoints[indices[i]]._x / mx;
            yc[i] = dpoints[indices[i]]._y / mx;
        }

        kdTree = new KDTree(0, xc, yc, indices);
    }

    int n;
    scanf("%d", &n);
    Integers result;
    Integers dindices;
    Integers temp;
    vector<bool> used;
    used.resize(m);
    queue<int> qu;

    for (int i = 0; i < n; ++i) {
        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        DPoint dq(x, y);

        dindices.clear();
        
        temp.clear();
        TBase mind = INF;
        for (int i = 0; i < 50; ++i) {
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
                for (Integers::const_iterator toItem = graph[index].begin(); toItem != graph[index].end(); ++toItem) {
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

        {
            float xn = x / mx;
            float yn = y / mx;

            TVector x4;
            x4.v = _mm_set1_ps(xn);
            TVector y4;
            y4.v = _mm_set1_ps(yn);

            float minDist = mind / mx / mx + 1e-6f;
            TVector minMax;
            minMax.v = _mm_set1_ps(minDist);
            float minMin = minDist;
            int limit = 10000;
            kdTree->Solve(xn, yn, x4, y4, &minDist, &minMax, &minMin, &dindices, &limit);
        }

        sort(dindices.begin(), dindices.end());

        result.clear();
        long double min = 1e15;
        static const long double LDEPS = 1e-10;
        int prevIndex = -1;
        for (int k = 0; k < dindices.size(); ++k) {
            int index = dindices[k];
            if (index > prevIndex) {
                long double d = dpoints[index].distance2(dq);
                if (d + LDEPS < min) {
                    result.clear();
                }
                if (d < min) {
                    min = d;
                }
                if (d < min + LDEPS) {
                    result.push_back(index);
                }
                prevIndex = index;
            }
        }

        Output(result);
    }

    return 0;
}