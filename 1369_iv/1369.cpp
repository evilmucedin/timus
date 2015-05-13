#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cmath>

#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <memory>
#include <iostream>

#define M_PI 3.14159265358979323846

using namespace std;

static const double INF = 1e9;

typedef vector<int> Integers;

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

void Output(const Integers& vct)
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

typedef double TBase;

struct Point {
    TBase _x;
    TBase _y;
    int _index;
    static const int ONSEGMENT = 0;
    static const int LEFT = 1;
    static const int RIGHT = 2;
    static const int INFRONTOFA = 3;
    static const int BEHINDB = 4;
    static const int ERROR = 5;

    Point()
        : _index(-1)
    {
    }

    Point(TBase x, TBase y, int index) {
        _x = x;
        _y = y;
        _index = index;
    }

    Point(const Point& p) {
        _x = p._x;
        _y = p._y;
        _index = p._index;
    }

    TBase x() const {
        return _x;
    }

    TBase y() const {
        return _y;
    }

    int index() const {
        return _index;
    }

    TBase distance2(const Point& p) const {
        return Sqr(p._x - _x) + Sqr(p._y - _y);
    }

    TBase distance2(TBase px, TBase py) const {
        return Sqr(px - _x) + Sqr(py - _y);
    }

    bool isLess(const Point& p) const {
        return (_x < p._x) || ((_x == p._x) && (_y < p._y));
    }

    bool operator<(const Point& p) const {
        return isLess(p);
    }

    bool isGreater(const Point& p) const {
        return (_x > p._x) || ((_x == p._x) && (_y > p._y));
    }

    bool equals(const Point& p) const {
        return (_x == p._x) && (_y == p._y);
    }

    bool operator==(const Point& rhs) const {
        return equals(rhs);
    }

    TBase distance(const Point& p) {
        TBase temp = Sqr(p.x() - _x) + Sqr(p.y() - _y);
        return sqrt(temp);
    }

    int pointLineTest(const Point& a, const Point& b) const {
        TBase dx = b._x - a._x;
        TBase dy = b._y - a._y;
        TBase res = dy * (_x - a._x) - dx * (_y - a._y);
        if (res < 0.0) {
            return 1;
        }
        if (res > 0.0) {
            return 2;
        }
        if (dx > 0.0) {
            if (_x < a._x) {
                return 3;
            }
            if (b._x < _x) {
                return 4;
            }
            return 0;
        }
        if (dx < 0.0) {
            if (_x > a._x) {
                return 3;
            }
            if (b._x > _x) {
                return 4;
            }
            return 0;
        }
        if (dy > 0.0) {
            if (_y < a._y) {
                return 3;
            }
            if (b._y < _y) {
                return 4;
            }
            return 0;
        }
        if (dy < 0.0) {
            if (_y > a._y) {
                return 3;
            }
            if (b._y > _y) {
                return 4;
            }
            return 0;
        }
        throw Exception("Error, pointLineTest with a=b");
    }

    bool areCollinear(const Point& a, const Point& b) {
        TBase dx = b._x - a._x;
        TBase dy = b._y - a._y;
        TBase res = dy * (_x - a._x) - dx * (_y - a._y);
        return res == 0.0;
    }

    bool areCollinearEps(const Point* a, const Point* b) {
        if (!a)
            return false;
        if (!b)
            return false;
        TBase dx = b->_x - a->_x;
        TBase dy = b->_y - a->_y;
        TBase res = dy * (_x - a->_x) - dx * (_y - a->_y);
        static const TBase EPS = 1e-4;
        return (res < EPS) && (res > -EPS);
    }

    Point circumcenter(const Point& a, const Point& b) {
        TBase u = ((a._x - b._x) * (a._x + b._x) + (a._y - b._y) * (a._y + b._y));
        TBase v = ((b._x - _x) * (b._x + _x) + (b._y - _y) * (b._y + _y));
        TBase den = ((a._x - b._x) * (b._y - _y) - (b._x - _x) * (a._y - b._y))*2.0;
        if (den == 0.0) {
            throw Exception("circumcenter, degenerate case");
        }
        return Point((u * (b._y - _y) - v * (a._y - b._y)) / den, (v * (a._x - b._x) - u * (b._x - _x)) / den, -1);
    }
};

typedef vector<Point> Points;

Points gPoints;

struct Circle {
    Point _c;
    TBase _r2;

    Circle() {
    }

    Circle(const Point& c, TBase r) {
        _c = c;
        _r2 = r;
    }

    Circle(const Circle& circ) {
        _c = circ._c;
        _r2 = circ._r2;
    }

    bool contains(const Point& p) const {
        return _c.distance2(p) < _r2;
    }

    Circle& operator=(const Circle& rhs) {
        _c = rhs._c;
        _r2 = rhs._r2;
        return *this;
    }
};

struct Triangle;
typedef shared_ptr<Triangle> PTriangle;

struct Triangle {
    int _a;
    int _b;
    int _c;
    PTriangle _abnext;
    PTriangle _bcnext;
    PTriangle _canext;
    bool _halfplane;
    bool _mark;
    Circle _circle;

    Triangle(int a, int b, int c)
        : _mark(false)
        , _halfplane(false)
    {
        _a = a;
        int res = gPoints[c].pointLineTest(gPoints[a], gPoints[b]);
        if ((res <= 1) || (res == 3) || (res == 4)) {
            _b = b;
            _c = c;
        }
        else {
            _b = c;
            _c = b;
        }
        cacheCircumCircle();
    }

    Triangle(int a, int b)
        : _mark(false)
    {
        _a = a;
        _b = b;
        _c = -1;
        _halfplane = true;
    }

    bool isHalfplane() {
        return _halfplane;
    }

    const Point& p1() const {
        return gPoints[_a];
    }

    const Point& p2() const {
        return gPoints[_b];
    }

    const Point& p3() const {
        return gPoints[_c];
    }

    PTriangle next_12() {
        return _abnext;
    }

    PTriangle next_23() {
        return _bcnext;
    }

    PTriangle next_31() {
        return _canext;
    }

    void switchNeighbors(PTriangle oldT, PTriangle newT) {
        if (_abnext == oldT) {
            _abnext = newT;
        }
        else if (_bcnext == oldT) {
            _bcnext = newT;
        }
        else if (_canext == oldT) {
            _canext = newT;
        }
        else {
            throw Exception("Error, switchneighbors can't find Old.");
        }
    }

    PTriangle neighbor(int p) const {
        if (_a == p) {
            return _canext;
        }
        if (_b == p) {
            return _abnext;
        }
        if (_c == p) {
            return _bcnext;
        }
        throw Exception("Error, neighbors can't find p: ");
    }

    Circle circumCircleRaw() const {
        TBase u = ((p1()._x - p2()._x) * (p1()._x + p2()._x) + (p1()._y - p2()._y) * (p1()._y + p2()._y));
        TBase v = ((p2()._x - p3()._x) * (p2()._x + p3()._x) + (p2()._y - p3()._y) * (p2()._y + p3()._y));
        TBase den = ((p1()._x - p2()._x) * (p2()._y - p3()._y) - (p2()._x - p3()._x) * (p1()._y - p2()._y))*2.0;
        if (den == 0.0) {
            return Circle(p1(), INF);
        } else {
            Point cen((u * (p2()._y - p3()._y) - v * (p1()._y - p2()._y)) / den, (v * (p1()._x - p2()._x) - u * (p2()._x - p3()._x)) / den, -1);
            return Circle(cen, cen.distance2(p1()));
        }
    }

    void cacheCircumCircle() {
        _circle = circumCircleRaw();
    }

    bool circumcircleContains(const Point& p) const {
        return _circle.contains(p);
    }

    bool contains(const Point* p) const {
        bool ans = false;
        if ((_halfplane || p == nullptr)) {
            return false;
        }
        if (p1().equals(*p) || p2().equals(*p) || p3().equals(*p)) {
            return true;
        }
        int a12 = p->pointLineTest(p1(), p2());
        int a23 = p->pointLineTest(p2(), p3());
        int a31 = p->pointLineTest(p3(), p1());
        if (((a12 == 1) && (a23 == 1) && (a31 == 1)) || ((a12 == 2) && (a23 == 2) && (a31 == 2)) || (a12 == 0) || (a23 == 0) || (a31 == 0)) {
            ans = true;
        }
        return ans;
    }
};

typedef vector<Point> Points;
typedef vector<DPoint> DPoints;

struct DelaunayTriangulation {
    int _firstP;
    int _lastP;
    bool _allCollinear;
    PTriangle _firstT;
    PTriangle _lastT;
    PTriangle _currT;
    PTriangle _startTriangle;
    PTriangle _startTriangleHull;
    int _nPoints;
    vector<PTriangle> _triangles;
    int _modCount;
    int _trianglesModCount;

    DelaunayTriangulation(const Integers& points)
        : _nPoints(0)
        , _modCount(0)
        , _trianglesModCount(0)
    {
        _modCount = 0;
        _trianglesModCount = 0;
        _allCollinear = true;
        for (int i = 0; i < points.size(); ++i) {
            insertPoint(points[i]);
        }
    }

    int trianglesSize() {
        initTriangles();
        return _triangles.size();
    }

    int getModificationsCounter() const {
        return _modCount;
    }

    void insertPoint(int p) {
        ++_modCount;
        PTriangle t = insertPointSimple(p);
        if (t == nullptr) {
            return;
        }
        PTriangle tt = t;
        _currT = t;
        do {
            flip(tt);
            tt = tt->_canext;
        } while ((tt != t) && (!tt->_halfplane));
    }

    PTriangle insertPointSimple(int p) {
        ++_nPoints;
        if (!_allCollinear) {
            PTriangle t = find(_startTriangle, &(gPoints[p]));
            if (t->_halfplane) {
                _startTriangle = extendOutside(t, p);
            } else {
                _startTriangle = extendInside(t, p);
            }
            return _startTriangle;
        }
        if (_nPoints == 1) {
            _firstP = p;
            return nullptr;
        }
        if (_nPoints == 2) {
            startTriangulation(_firstP, p);
            return nullptr;
        }
        switch (gPoints[p].pointLineTest(gPoints[_firstP], gPoints[_lastP])) {
        case 1:
            _startTriangle = extendOutside(_firstT->_abnext, p);
            _allCollinear = false;
            break;
        case 2:
            _startTriangle = extendOutside(_firstT, p);
            _allCollinear = false;
            break;
        case 0:
            insertCollinear(p, 0);
            break;
        case 3:
            insertCollinear(p, 3);
            break;
        case 4:
            insertCollinear(p, 4);
        }
        return nullptr;
    }

    void insertCollinear(int p, int res) {
        switch (res) {
        case 3: {
                    PTriangle t(make_shared<Triangle>(_firstP, p));
                    PTriangle tp(make_shared<Triangle>(p, _firstP));
                    t->_abnext = tp;
                    tp->_abnext = t;
                    t->_bcnext = tp;
                    tp->_canext = t;
                    t->_canext = _firstT;
                    _firstT->_bcnext = t;
                    tp->_bcnext = _firstT->_abnext;
                    _firstT->_abnext->_canext = tp;
                    _firstT = t;
                    _firstP = p;
        }
            break;
        case 4: {
                    PTriangle t(make_shared<Triangle>(p, _lastP));
                    PTriangle tp(make_shared<Triangle>(_lastP, p));
                    t->_abnext = tp;
                    tp->_abnext = t;
                    t->_bcnext = _lastT;
                    _lastT->_canext = t;
                    t->_canext = tp;
                    tp->_bcnext = t;
                    tp->_canext = _lastT->_abnext;
                    _lastT->_abnext->_bcnext = tp;
                    _lastT = t;
                    _lastP = p;
        }
            break;
        case 0: {
                    PTriangle u = _firstT;
                    while (gPoints[p].isGreater(u->p1())) {
                        u = u->_canext;
                    }
                    PTriangle t(make_shared<Triangle>(p, u->_b));
                    PTriangle tp(make_shared<Triangle>(u->_b, p));
                    u->_b = p;
                    u->_abnext->_a = p;
                    t->_abnext = tp;
                    tp->_abnext = t;
                    t->_bcnext = u->_bcnext;
                    u->_bcnext->_canext = t;
                    t->_canext = u;
                    u->_bcnext = t;
                    tp->_canext = u->_abnext->_canext;
                    u->_abnext->_canext->_bcnext = tp;
                    tp->_bcnext = u->_abnext;
                    u->_abnext->_canext = tp;
                    if (_firstT == u) {
                        _firstT = t;
                    }
        }
            break;
        }
    }

    void startTriangulation(int p1, int p2) {
        int ps;
        int pb;
        if (gPoints[p1].isLess(gPoints[p2])) {
            ps = p1;
            pb = p2;
        } else {
            ps = p2;
            pb = p1;
        }
        _firstT = make_shared<Triangle>(pb, ps);
        _lastT = _firstT;
        PTriangle t(make_shared<Triangle>(ps, pb));
        _firstT->_abnext = t;
        t->_abnext = _firstT;
        _firstT->_bcnext = t;
        t->_canext = _firstT;
        _firstT->_canext = t;
        t->_bcnext = _firstT;
        _firstP = _firstT->_b;
        _lastP = _lastT->_a;
        _startTriangleHull = _firstT;
    }

    PTriangle extendInside(PTriangle t, int p) {
        PTriangle h1 = treatDegeneracyInside(t, p);
        if (h1 != nullptr) {
            return h1;
        }
        h1 = PTriangle(make_shared<Triangle>(t->_c, t->_a, p));
        PTriangle h2(make_shared<Triangle>(t->_b, t->_c, p));
        t->_c = p;
        t->cacheCircumCircle();
        h1->_abnext = t->_canext;
        h1->_bcnext = t;
        h1->_canext = h2;
        h2->_abnext = t->_bcnext;
        h2->_bcnext = h1;
        h2->_canext = t;
        h1->_abnext->switchNeighbors(t, h1);
        h2->_abnext->switchNeighbors(t, h2);
        t->_bcnext = h2;
        t->_canext = h1;
        return t;
    }

    PTriangle treatDegeneracyInside(PTriangle t, int p) {
        if ((t->_abnext->_halfplane) && (gPoints[p].pointLineTest(t->p2(), t->p1()) == 0)) {
            return extendOutside(t->_abnext, p);
        }
        if ((t->_bcnext->_halfplane) && (gPoints[p].pointLineTest(t->p3(), t->p2()) == 0)) {
            return extendOutside(t->_bcnext, p);
        }
        if ((t->_canext->_halfplane) && (gPoints[p].pointLineTest(t->p1(), t->p3()) == 0)) {
            return extendOutside(t->_canext, p);
        }
        return nullptr;
    }

    PTriangle extendOutside(PTriangle t, int p) {
        if (gPoints[p].pointLineTest(t->p1(), t->p2()) == 0) {
            PTriangle dg(make_shared<Triangle>(t->_a, t->_b, p));
            PTriangle hp(make_shared<Triangle>(p, t->_b));
            t->_b = p;
            t->cacheCircumCircle();
            dg->_abnext = t->_abnext;
            dg->_abnext->switchNeighbors(t, dg);
            dg->_bcnext = hp;
            hp->_abnext = dg;
            dg->_canext = t;
            t->_abnext = dg;
            hp->_bcnext = t->_bcnext;
            hp->_bcnext->_canext = hp;
            hp->_canext = t;
            t->_bcnext = hp;
            return dg;
        }
        PTriangle ccT = extendCounterClock(t, p);
        PTriangle cT = extendClock(t, p);
        ccT->_bcnext = cT;
        cT->_canext = ccT;
        _startTriangleHull = cT;
        return cT->_abnext;
    }

    PTriangle extendCounterClock(PTriangle t, int p) {
        t->_halfplane = false;
        t->_c = p;
        t->cacheCircumCircle();

        PTriangle tca = t->_canext;
        if (gPoints[p].pointLineTest(tca->p1(), tca->p2()) >= 2) {
            PTriangle nT(make_shared<Triangle>(t->_a, p));
            nT->_abnext = t;
            t->_canext = nT;
            nT->_canext = tca;
            tca->_bcnext = nT;
            return nT;
        }
        return extendCounterClock(tca, p);
    }

    PTriangle extendClock(PTriangle t, int p) {
        t->_halfplane = false;
        t->_c = p;
        t->cacheCircumCircle();

        PTriangle tbc = t->_bcnext;
        if (gPoints[p].pointLineTest(tbc->p1(), tbc->p2()) >= 2) {
            PTriangle nT(make_shared<Triangle>(p, t->_b));
            nT->_abnext = t;
            t->_bcnext = nT;
            nT->_bcnext = tbc;
            tbc->_canext = nT;
            return nT;
        }
        return extendClock(tbc, p);
    }

    void flip(PTriangle t) {
        PTriangle u = t->_abnext;
        if ((u->_halfplane) || (!u->circumcircleContains(t->p3()))) {
            return;
        }
        PTriangle v;
        if (t->_a == u->_a) {
            v = PTriangle(make_shared<Triangle>(u->_b, t->_b, t->_c));
            v->_abnext = u->_bcnext;
            t->_abnext = u->_abnext;
        }
        else if (t->_a == u->_b) {
            v = PTriangle(make_shared<Triangle>(u->_c, t->_b, t->_c));
            v->_abnext = u->_canext;
            t->_abnext = u->_bcnext;
        }
        else if (t->_a == u->_c) {
            v = PTriangle(make_shared<Triangle>(u->_a, t->_b, t->_c));
            v->_abnext = u->_abnext;
            t->_abnext = u->_canext;
        }
        else {
            throw Exception("Error in flip.");
        }
        v->_bcnext = t->_bcnext;
        v->_abnext->switchNeighbors(u, v);
        v->_bcnext->switchNeighbors(t, v);
        t->_bcnext = v;
        v->_canext = t;
        t->_b = v->_a;
        t->cacheCircumCircle();
        t->_abnext->switchNeighbors(u, t);

        _currT = v;
        flip(t);
        flip(v);
    }

    PTriangle find(const Point* p) {
        return find(_startTriangle, p);
    }

    PTriangle find(const Point* p, PTriangle start) {
        if (start == nullptr) {
            start = _startTriangle;
        }
        return find(start, p);
    }

    PTriangle find(PTriangle curr, const Point* p) {
        if (p == nullptr) {
            return nullptr;
        }
        if (curr->_halfplane) {
            PTriangle next_t = findnext2(*p, curr);
            if ((next_t == nullptr) || (next_t->_halfplane)) {
                return curr;
            }
            curr = next_t;
        }
        for (;;) {
            PTriangle next_t = findnext1(*p, curr);
            if (next_t == nullptr) {
                return curr;
            }
            if (next_t->_halfplane) {
                return next_t;
            }
            curr = next_t;
        }
    }

    PTriangle findnext1(const Point& p, PTriangle v) {
        if ((p.pointLineTest(v->p1(), v->p2()) == 2) && (!v->_abnext->_halfplane)) {
            return v->_abnext;
        }
        if ((p.pointLineTest(v->p2(), v->p3()) == 2) && (!v->_bcnext->_halfplane)) {
            return v->_bcnext;
        }
        if ((p.pointLineTest(v->p3(), v->p1()) == 2) && (!v->_canext->_halfplane)) {
            return v->_canext;
        }
        if (p.pointLineTest(v->p1(), v->p2()) == 2) {
            return v->_abnext;
        }
        if (p.pointLineTest(v->p2(), v->p3()) == 2) {
            return v->_bcnext;
        }
        if (p.pointLineTest(v->p3(), v->p1()) == 2) {
            return v->_canext;
        }
        return nullptr;
    }

    PTriangle findnext2(const Point& p, PTriangle v) {
        if ((v->_abnext != nullptr) && (!v->_abnext->_halfplane)) {
            return v->_abnext;
        }
        if ((v->_bcnext != nullptr) && (!v->_bcnext->_halfplane)) {
            return v->_bcnext;
        }
        if ((v->_canext != nullptr) && (!v->_canext->_halfplane)) {
            return v->_canext;
        }
        return nullptr;
    }

    bool contains(const Point& p) {
        PTriangle tt = find(&p);
        return !tt->_halfplane;
    }

    PTriangle safeFind(const Point& q) {
        trianglesIterator();
        for (size_t i = 0; i < _triangles.size(); ++i) {
            if (_triangles[i]->contains(&q)) {
                return _triangles[i];
            }
        }
        return nullptr;
    }

    void trianglesIterator() {
        if (_nPoints > 2) {
            initTriangles();
        }
    }

    void initTriangles() {
        if (_modCount == _trianglesModCount) {
            return;
        }
        if (_nPoints > 2) {
            _trianglesModCount = _modCount;
            queue<PTriangle> front;
            _triangles.clear();
            front.push(_startTriangle);
            while (!front.empty()) {
                PTriangle t = front.front();
                front.pop();
                if (!t->_mark) {
                    t->_mark = true;
                    _triangles.push_back(t);
                    if ((t->_abnext != nullptr) && (!t->_abnext->_mark)) {
                        front.push(t->_abnext);
                    }
                    if ((t->_bcnext != nullptr) && (!t->_bcnext->_mark)) {
                        front.push(t->_bcnext);
                    }
                    if ((t->_canext != nullptr) && (!t->_canext->_mark)) {
                        front.push(t->_canext);
                    }
                }
            }
            for (int i = 0; i < _triangles.size(); ++i) {
                _triangles[i]->_mark = false;
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
    for (int i = 0; i < M; ++i)
    {
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

int main() {
#ifndef ONLINE_JUDGE
    // GenBig();
    // freopen("big.txt", "r", stdin);
    freopen("input.txt", "r", stdin);
#endif

    // cout << sizeof(Triangle) << endl;

    int m;
    scanf("%d", &m);
    Points points(m);
    DPoints dpoints(m);
    Integers indices(m);
    for (int i = 0; i < m; ++i) {
        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        points[i] = Point(x, y, i);
        dpoints[i] = DPoint(x, y);
        indices[i] = i;
    }

    for (int i = 0; i < m; ++i)
    {
        swap(indices[i], indices[i + (rand() % (m - i))]);
    }
    if (m > 40000)
        indices.erase(indices.begin() + 40000, indices.end());

    gPoints = points;
    DelaunayTriangulation dt(indices);

    vector<vector<int>> graph(m, vector<int>());
    dt.trianglesIterator();
    for (size_t i = 0; i < dt._triangles.size(); ++i)
    {
        PTriangle t = dt._triangles[i];
        if (t->_a >= 0 && t->_b >= 0) {
            graph[t->p1().index()].push_back(t->p2().index());
            graph[t->p2().index()].push_back(t->p1().index());
        }

        if (t->_a >= 0 && t->_c >= 0) {
            graph[t->p1().index()].push_back(t->p3().index());
            graph[t->p3().index()].push_back(t->p1().index());
        }

        if (t->_b >= 0 && t->_c >= 0) {
            graph[t->p2().index()].push_back(t->p3().index());
            graph[t->p3().index()].push_back(t->p2().index());
        }
    }

    int n;
    scanf("%d", &n);
    vector<bool> used(m);
    Integers result;
    queue<int> qu;

    if (m > 10000)
        return 0;
    
    for (int i = 0; i < n; ++i) {        
        indices.clear();

        long double x;
        long double y;
        scanf("%Lf%Lf", &x, &y);
        Point q(x, y, -1);
        DPoint dq(x, y);

        bool fallback = false;
        if (m > 5) {
            PTriangle t = dt.find(&q);
            if (nullptr != t) {
                if (t->_a >= 0) {
                    indices.push_back(t->p1().index());
                }
                if (t->_b >= 0) {
                    indices.push_back(t->p2().index());
                }
                if (t->_c >= 0) {
                    indices.push_back(t->p3().index());
                }

                for (int k = 0; k < indices.size(); ++k) {
                    qu.push(indices[k]);
                    used[indices[k]] = true;
                }
                TBase min = INF;
                while (!qu.empty()) {
                    int index = qu.front();
                    qu.pop();
                    double dist = points[index].distance2(q);
                    if (dist < min + 1e-9) {
                        if (dist < min) {
                            min = dist;
                        }
                        for (size_t k = 0; k < graph[index].size(); ++k) {
                            int v = graph[index][k];
                            if (!used[v]) {
                                indices.push_back(v);
                                used[v] = true;
                                qu.push(v);
                            }
                        }
                    }
                }
                for (int k = 0; k < indices.size(); ++k) {
                    used[indices[k]] = false;
                }
            }
            else {
                fallback = true;
            }
        }
        else {
            fallback = true;
        }
        if (fallback) {
            indices.resize(m);
            for (int j = 0; j < m; ++j) {
                indices[j] = j;
            }
        }
        else {
            sort(indices.begin(), indices.end());
        }

        // BigDecimal min = new BigDecimal(1e12);
        result.clear();
        long double min = 1e12;
        int prevIndex = -1;
        for (int k = 0; k < indices.size(); ++k) {
            int index = indices[k];
            if (index > prevIndex) {
                long double d = points[index].distance2(q);
                // int cmp = d.compareTo(min);
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