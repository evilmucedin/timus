#define _CRT_SECURE_NO_WARNINGS

#ifndef _MSC_VER
#   pragma GCC target("sse4.2")
#   pragma GCC optimize("O3")
#endif

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

static const long double INF = 1e30;
static const long double EPS = 5e-12;

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

typedef long double TBase;

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

    TBase Distance2(const Point& p) const {
        return Sqr(_x - p._x) + Sqr(_y - p._y);
    }

    bool operator==(const Point& rhs) {
        return Distance2(rhs) < EPS;
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
    __m128d v;
    double f[2];
};

typedef vector<long double> TFloats;

static const float FINF = 1e8f;
static const double DINF = 1e10;

static void ConvertVector(const TFloats& floats, TVector** result, int* len)
{
    const int m = static_cast<int>(floats.size());
    *len = m / 2;
    if (m % 2)
    {
        ++(*len);
    }
    *result = (TVector*)_mm_malloc(sizeof(TVector)*(*len), 32);
    for (int i = 0; i < *len; ++i)
    {
        double values[] = { DINF, DINF };
        for (int j = 0; j < 2; ++j)
        {
            if (2 * i + j < m)
            {
                values[j] = floats[2 * i + j];
            }
        }
        (*result)[i].v = _mm_load_pd(values);
    }
}

static long double Median(const TFloats& floats) {
    TFloats temp(floats);
    size_t n = temp.size() / 2;
    nth_element(temp.begin(), temp.begin() + n, temp.end());
    return temp[n];
}

struct KDTree
{
    int _len;
    int _depth;
    TVector* _x;
    TVector* _y;
    bool _isLeaf;
    KDTree* _left;
    KDTree* _right;
    Integers _indices;
    long double _minX;
    long double _maxX;
    long double _minY;
    long double _maxY;
    long double _separator;

    KDTree(int depth, const TFloats& x, const TFloats& y, const Integers& indices)
        : _left(nullptr)
        , _right(nullptr)
    {
        _minX = FINF;
        _maxX = -FINF;
        _minY = FINF;
        _maxY = -FINF;
        for (size_t i = 0; i < x.size(); ++i) {
            _minX = Min(_minX, x[i]);
            _maxX = Max(_maxX, x[i]);
            _minY = Min(_minY, y[i]);
            _maxY = Max(_maxY, y[i]);
        }
        _depth = depth;

        if (12 == depth) {
            ConvertVector(x, &_x, &_len);
            ConvertVector(y, &_y, &_len);
            _indices = indices;
            while (_indices.size() < 2 * _len) {
                _indices.push_back(-1);
            }
            _isLeaf = true;
        } else {
            _separator = (depth & 1) ? Median(x) : Median(y);

            TFloats xLeft;
            TFloats yLeft;
            Integers indicesLeft;
            TFloats xRight;
            TFloats yRight;
            Integers indicesRight;

            for (size_t i = 0; i < x.size(); ++i) {
                bool toLeft = (depth & 1) ? (x[i] < _separator) : (y[i] < _separator);
                if (toLeft) {
                    xLeft.push_back(x[i]);
                    yLeft.push_back(y[i]);
                    indicesLeft.push_back(indices[i]);
                } else {
                    xRight.push_back(x[i]);
                    yRight.push_back(y[i]);
                    indicesRight.push_back(indices[i]);
                }
            }

            if (xLeft.size()) {
                _left = new KDTree(depth + 1, xLeft, yLeft, indicesLeft);
            }
            if (xRight.size()) {
                _right = new KDTree(depth + 1, xRight, yRight, indicesRight);
            }

            _isLeaf = false;
        }
    }

#ifndef _MSC_VER
    __attribute__((force_align_arg_pointer))
#endif
        void Solve(long double x, long double y, const TVector& x2, const TVector& y2, long double* minDist, long double* minDistEps, TVector* minMax, Integers* result, int* limit) const
    {
        if (*limit < 0) {
            return;
        }

        long double d2;
        if (x < _minX) {
            if (y < _minY) {
                d2 = Sqr(_minX - x) + Sqr(_minY - y);
            } else if (y > _maxY) {
                d2 = Sqr(_minX - x) + Sqr(y - _maxY);
            } else {
                d2 = Sqr(_minX - x);
            }
        } else if (x > _maxX) {
            if (y < _minY) {
                d2 = Sqr(x - _maxX) + Sqr(_minY - y);
            } else if (y > _maxY) {
                d2 = Sqr(x - _maxX) + Sqr(y - _maxY);
            } else {
                d2 = Sqr(x - _maxX);
            }
        } else {
            if (y < _minY) {
                d2 = Sqr(_minY - y);
            } else if (y > _maxY) {
                d2 = Sqr(y - _maxY);
            } else {
                d2 = 0;
            }
        }

        if (d2 > *minDist + EPS) {
            return;
        }

        if (_isLeaf) {
            *limit -= _len;
            for (int j = 0; j < _len; ++j) {
                TVector dx2;
                dx2.v = _mm_sub_pd(_x[j].v, x2.v);
                dx2.v = _mm_mul_pd(dx2.v, dx2.v);
                TVector dy2;
                dy2.v = _mm_sub_pd(_y[j].v, y2.v);
                dy2.v = _mm_mul_pd(dy2.v, dy2.v);
                dx2.v = _mm_add_pd(dx2.v, dy2.v);

                TVector cmpMax;
                cmpMax.v = _mm_cmple_pd(dx2.v, minMax->v);
                for (int k = 0; k < 2; ++k) {
                    if (cmpMax.f[k]) {
                        if (dx2.f[k] <= *minDistEps) {
                            if (dx2.f[k] < *minDist) {
                                if (dx2.f[k] + EPS < *minDist) {
                                    result->clear();
                                }
                                *minDist = dx2.f[k];
                                *minDistEps = *minDist + EPS;
                                minMax->v = _mm_set1_pd(*minDistEps);
                            }
                            result->push_back(_indices[2 * j + k]);
                        }
                    }
                }
            }
        } else {
            bool toLeft = (_depth & 1) ? (x < _separator) : (y < _separator);
            if (toLeft) {
                if (_left) {
                    _left->Solve(x, y, x2, y2, minDist, minDistEps, minMax, result, limit);
                }
                if (_right) {
                    _right->Solve(x, y, x2, y2, minDist, minDistEps, minMax, result, limit);
                }
            } else {
                if (_right) {
                    _right->Solve(x, y, x2, y2, minDist, minDistEps, minMax, result, limit);
                }
                if (_left) {
                    _left->Solve(x, y, x2, y2, minDist, minDistEps, minMax, result, limit);
                }
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

#ifndef _MSC_VER
#   define GET_ITEM2(var, index) (var)[(index)]
#else
#   define GET_ITEM2(var, index) (var).m128d_f64[(index)]
#endif

struct QDouble {
    long double _int;
    long double _float;

    QDouble() {

    }

    static const long double EPS;
    static const long double NORM;

    QDouble(long double ld) {
        ld *= NORM;
        _int = floor(ld);
        _float = ld - _int;
        Normalize();
    }

    void Normalize() {
        long double diff2 = _int - floor(_int);
        _int -= diff2;
        _float += diff2;

        long double diff = floor(_float);
        _float -= diff;
        _int += diff;

        if (_float < 0.0) {
            _float += 1.0;
            _int -= 1.0;
        }
        if (_float > 1.0) {
            _float -= 1.0;
            _int += 1.0;
        }
    }

    long double ToDouble() const {
        return (_int + _float)/NORM;
    }

    bool operator<(const QDouble& rhs) const {
        if (_int < rhs._int) {
            return true;
        }
        else if (_int == rhs._int) {
            return _float < rhs._float;
        }
        return false;
    }

    bool operator<=(const QDouble& rhs) const {
        if (_int == rhs._int) {
            return _float <= rhs._float;
        }
        else if (_int < rhs._int) {
            return true;
        }
        return false;
    }

    void DeNorm() {
        _int /= NORM;
        _float /= NORM;
    }
};

const long double QDouble::EPS = 1e-12;
const long double QDouble::NORM = 100.0;

void Minus(const QDouble& q1, const QDouble& q2, QDouble* result) {
    result->_int = q1._int - q2._int;
    result->_float = q1._float - q2._float;
    result->Normalize();
}

void Add(const QDouble& q1, const QDouble& q2, QDouble* result) {
    result->_int = q1._int + q2._int;
    result->_float = q1._float + q2._float;
    result->Normalize();
}

void Sqr(const QDouble& q, QDouble* result) {
    result->_int = q._int*q._int;
    result->_float = 2.0*q._float*q._int + q._float*q._float;
    result->DeNorm();
    result->Normalize();
}

void Mul(const QDouble& a, const QDouble& b, QDouble* result) {
    result->_int = a._int*b._int;
    result->_float = a._float*b._int + b._float*a._int + a._float*b._float;
    result->DeNorm();
    result->Normalize();
}

void ReadQDouble(const char* s, QDouble* result) {
    QDouble power10(1.0);
    QDouble now(0.0);
    static const QDouble ten(10.0);
    static const QDouble tenth(0.1);
    bool minus = false;
    bool inInt = true;
    while (*s) {
        if (*s == '-') {
            minus = true;
        }
        else if (*s == '.') {
            inInt = false;
        }
        else {
            int digit = *s - '0';
            if (inInt) {
                QDouble temp1;
                Mul(now, ten, &temp1);
                QDouble temp2;
                Add(temp1, QDouble(digit), &temp2);
                now = temp2;
            } else {
                QDouble temp1;
                Mul(power10, tenth, &temp1);
                power10 = temp1;

                QDouble temp2;
                Mul(power10, QDouble(digit), &temp2);

                QDouble temp3;
                Add(now, temp2, &temp3);

                now = temp3;
            }
        }
        ++s;
    }
    *result = now;
    if (minus) {
        result->_float = -result->_float;
        result->_int = -result->_int;
        result->Normalize();
    }
}

struct QPoint {
    QDouble _x;
    QDouble _y;

    QPoint() {

    }

    QPoint(long double x, long double y)
        : _x(x)
        , _y(y)
    {

    }

    QDouble distance2(const QPoint& p) const {
        QDouble dx;
        Minus(_x, p._x, &dx);
        QDouble dx2;
        Sqr(dx, &dx2);

        QDouble dy;
        Minus(_y, p._y, &dy);
        QDouble dy2;
        Sqr(dy, &dy2);

        QDouble sum;
        Add(dx2, dy2, &sum);

        return sum;
    }
};
typedef vector<QPoint> QPoints;

Point CircumCenter(const Point& a, const Point& b, const Point& c) {
    TBase u = ((a._x - b._x) * (a._x + b._x) + (a._y - b._y) * (a._y + b._y));
    TBase v = ((b._x - c._x) * (b._x + c._x) + (b._y - c._y) * (b._y + c._y));
    TBase den = ((a._x - b._x) * (b._y - c._y) - (b._x - c._x) * (a._y - b._y))*2.0;
    return Point((u * (b._y - c._y) - v * (a._y - b._y)) / den, (v * (a._x - b._x) - u * (b._x - c._x)) / den, -1);
}

struct Circle {
    Point _center;
    TBase _r2;
    Integers _indices;

    void MinD(const Point& q, TBase* mind, Integers* result) const {
        if (_r2 < *mind + EPS) {
            if (q.Distance2(_center) < EPS) {
                if (_r2 < *mind) {
                    if (_r2 + EPS < *mind) {
                        result->clear();
                    }
                    *mind = _r2;
                }
                result->insert(result->end(), _indices.begin(), _indices.end());
            }
        }
    }

    bool operator!=(const Circle& rhs) const {
        if (_center.Distance2(rhs._center) > EPS) {
            return true;
        }
        if (abs(_r2 - rhs._r2) > EPS) {
            return true;
        }
        return false;
    }
};

int main() {
#ifndef ONLINE_JUDGE
    // GenBig();
    // freopen("big.txt", "r", stdin);
    freopen("input.txt", "r", stdin);
#endif

    int m;
    scanf("%d", &m);
    DPoints dpoints(m);
    // QPoints qpoints(m);
    __m128d* xyd = (__m128d*)_mm_malloc(sizeof(__m128d)*(m + 8), 32);
    long double mx = 1e-12;
    char sx[1000];
    char sy[1000];
    for (int i = 0; i < m; ++i) {
        scanf("%s%s", sx, sy);
        // ReadQDouble(sx, &qpoints[i]._x);
        // ReadQDouble(sy, &qpoints[i]._y);

        long double x;
        sscanf(sx, "%Lf", &x);
        long double y;
        sscanf(sy, "%Lf", &y);

        dpoints[i] = DPoint(x, y);
        mx = Max(mx, abs(x));
        mx = Max(mx, abs(y));
        xyd[i] = _mm_set_pd(x, y);
    }
    mx = 1.0;

    Integers indices(m);
    for (int i = 0; i < m; ++i) {
        indices[i] = i;
    }
    for (int i = 0; i < m; ++i) {
        swap(indices[i], indices[i + (rand() % (m - i))]);
    }

    vector<Integers> graph;
    graph.resize(m);

    static const size_t NCENTERS = 4;
    Points centers;
    vector<PointProjections> angles;

    vector<Circle> circles;

    Points points(m);
    {
        for (int i = 0; i < m; ++i) {
            points[i] = Point(dpoints[indices[i]]._x, dpoints[indices[i]]._y, indices[i]);
        }

        centers.resize(NCENTERS);
        angles.resize(NCENTERS);
        for (int it = 0; it < NCENTERS; ++it) {
            Point cc = CircumCenter(points[rand() % m], points[rand() % m], points[rand() % m]);
            int iit = 0;
            while (iit < 100 && find(centers.begin(), centers.end(), cc) != centers.end() ) {
                cc = CircumCenter(points[rand() % m], points[rand() % m], points[rand() % m]);
                ++iit;
            }
            centers[it] = cc;
            angles[it].resize(m);

            for (int i = 0; i < m; ++i) {
                angles[it][i]._p = &points[i];
                angles[it][i]._projection = points[i].Distance2(centers[it]);
            }
            sort(angles[it].begin(), angles[it].end());

            int index = 0;
            while (index < angles[it].size()) {
                int endindex = index + 50;
                if (endindex < angles[it].size() && angles[it][endindex]._projection <= angles[it][index]._projection + EPS) {
                    while (endindex < angles[it].size() && angles[it][endindex]._projection <= angles[it][index]._projection + EPS) {
                        ++endindex;
                    }
                    Circle c;
                    c._center = centers[it];
                    c._r2 = angles[it][index]._projection;
                    int i = 0;
                    while (i < circles.size() && c != circles[i]) {
                        ++i;
                    }
                    if (i == circles.size()) {
                        c._indices.reserve(endindex - index);
                        for (int j = index; j < endindex; ++j) {
                            c._indices.push_back(angles[it][j]._p->_index);
                        }
                        circles.push_back(c);
                    }
                    index = endindex;
                } else {
                    ++index;
                }
            }

            for (int i = 0; i < m; ++i) {
                angles[it][i]._p = &points[i];
                angles[it][i]._projection = Angle(points[i], centers[it]);
            }
            sort(angles[it].begin(), angles[it].end());
        }

        PointProjections projections(m);

        static const int N_IT = 3;
        static const int D = 1;

        {
            static const int N_DIT = 10;
            PointProjections dists(N_DIT);

            for (int i = 0; i < m; ++i) {
                for (int it = 0; it < N_DIT; ++it) {
                    int index = rand() % m;
                    dists[it]._p = &points[index];
                    dists[it]._projection = points[index].Distance2(points[i]);
                }
                sort(dists.begin(), dists.end());
                int it = 0;
                while (it < dists.size() && (dists[it]._projection < EPS)) {
                    ++it;
                }
                if (it < dists.size()) {
                    int index1 = points[i]._index;
                    while (it < D && it < dists.size()) {
                        int index2 = dists[it]._p->_index;
                        graph[index1].push_back(index2);
                        graph[index2].push_back(index1);
                        ++it;
                    }
                }
            }
        }

        static const Vector vcts[2] = { Vector(1.0, 0.0), Vector(0.0, 1.0) };

        for (int it = -2; it < N_IT; ++it) {
            Vector vRand(Rand(), Rand());
            if (it == -2)
                vRand = vcts[0];
            else if (it == -1)
                vRand = vcts[1];
            vRand.Normalize();
            for (int i = 0; i < m; ++i) {
                projections[i]._p = &points[i];
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

        for (int it = 0; it < N_IT; ++it) {
            int index = rand() % m;
            for (int i = 0; i < m; ++i) {
                projections[i]._p = &points[i];
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
        TFloats xc(m);
        TFloats yc(m);
        for (int i = 0; i < m; ++i) {
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

    Integers components;
    for (int i = 0; i < m; ++i) {
        if (!used[i]) {
            components.push_back(i);

            qu.push(i);
            while (!qu.empty()) {
                int index = qu.front();
                qu.pop();
                for (Integers::const_iterator toItem = graph[index].begin(); toItem != graph[index].end(); ++toItem) {
                    int v = *toItem;
                    if (!used[v]) {
                        used[v] = true;
                        qu.push(v);
                    }
                }
            }
        }
    }
    for (int i = 0; i < m; ++i) {
        used[i] = false;
    }

    for (int i = 0; i < n; ++i) {
        QPoint qq;
        scanf("%s%s", sx, sy);
        ReadQDouble(sx, &qq._x);
        ReadQDouble(sy, &qq._y);

        long double x;
        sscanf(sx, "%Lf", &x);
        long double y;
        sscanf(sy, "%Lf", &y);

        DPoint dq(x, y);
        Point pq(x, y, -1);
        __m128d xy = _mm_set_pd(x, y);

        dindices.clear();

        TBase mind = INF;

        temp.clear();
        for (int j = 0; j < components.size(); ++j) {
            int index = components[j];
            qu.push(index);
            used[index] = true;
            temp.push_back(index);
        }
        for (int j = 0; j < 200; ++j) {
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

            __m128d xyd2 = _mm_sub_pd(xyd[index], xy);
            xyd2 = _mm_mul_pd(xyd2, xyd2);
            long double dist = GET_ITEM2(xyd2, 0) + GET_ITEM2(xyd2, 1);

            if (dist <= mind + EPS) {
                if (dist < mind) {
                    if (dist + EPS < mind) {
                        dindices.clear();
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

        for (int k = 0; k < angles.size(); ++k) {
            PointProjection qpp;
            qpp._p = &pq;
            qpp._projection = Angle(pq, centers[k]);
            vector<PointProjection>::const_iterator toAngle = lower_bound(angles[k].begin(), angles[k].end(), qpp);
            const PointProjection* bgn = &angles[k][0];
            for (int delta = -50; delta <= 50; ++delta) {
                const PointProjection* toCand = &(*toAngle) + delta;
                while (toCand < bgn)
                    toCand += angles[k].size();
                while (toCand >= bgn + angles[k].size())
                    toCand -= angles[k].size();

                const int index = toCand->_p->_index;
                __m128d xyd2 = _mm_sub_pd(xyd[index], xy);
                xyd2 = _mm_mul_pd(xyd2, xyd2);
                long double d = GET_ITEM2(xyd2, 0) + GET_ITEM2(xyd2, 1);

                if (d <= mind + EPS) {
                    if (d < mind) {
                        if (d + EPS < mind) {
                            dindices.clear();
                        }
                        mind = d;
                    }
                    dindices.push_back(index);
                }
            }
        }

        for (int k = 0; k < circles.size(); ++k) {
            circles[k].MinD(pq, &mind, &dindices);
        }

        long double xn = x / mx;
        long double yn = y / mx;

        TVector x2;
        x2.v = _mm_set1_pd(xn);
        TVector y2;
        y2.v = _mm_set1_pd(yn);

        long double minDist = mind / mx / mx;
        long double minDistEps = minDist + EPS;
        TVector minMax;
        minMax.v = _mm_set1_pd(minDistEps);
        int limit = 6000;
        kdTree->Solve(xn, yn, x2, y2, &minDist, &minDistEps, &minMax, &dindices, &limit);

        sort(dindices.begin(), dindices.end());
        dindices.erase(unique(dindices.begin(), dindices.end()), dindices.end());

        /*
        result.clear();
        long double min = 1e100;
        long double minPlus = 1e100;
        static const long double LDEPS(5e-12);
        int prevIndex = -1;
        for (int k = 0; k < dindices.size(); ++k) {
            int index = dindices[k];
            if (index > prevIndex) {
                // __m128d xyd2 = _mm_sub_pd(xyd[index], xy);
                // xyd2 = _mm_mul_pd(xyd2, xyd2);
                // long double d = GET_ITEM2(xyd2, 0) + GET_ITEM2(xyd2, 1);
                long double d = dpoints[index].distance2(dq);
                if (d + LDEPS < min) {
                    if (result.size()) {
                        throw Exception("error");
                    }
                    result.clear();
                }
                if (d < min) {
                    min = d;
                }
                if (d <= min + LDEPS) {
                    result.push_back(index);
                } else {
                    throw Exception("error");
                }
                prevIndex = index;
            }
        }
        */

        /*
        QDouble min = 1e100;
        QDouble minPlus = 1e100;
        static const QDouble LDEPS(1e-12);
        int prevIndex = -1;
        for (int k = 0; k < dindices.size(); ++k) {
        int index = dindices[k];
        if (index > prevIndex) {
        QDouble d = qq.distance2(qpoints[index]);
        QDouble temp;
        QDouble dPlus;
        Add(d, LDEPS, &dPlus);
        if (dPlus < min) {
        result.clear();
        }
        if (d < min) {
        min = d;
        Add(min, LDEPS, &minPlus);
        }
        if (d <= minPlus) {
        result.push_back(index);
        }
        prevIndex = index;
        }
        }
        */

        Output(dindices);
    }

    return 0;
}