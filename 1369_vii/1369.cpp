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

struct Edge;

struct Point {
    double x, y;
    int index;

    Point()
    {
        x = 0.0;
        y = 0.0;
        index = -1;
    }

    std::vector<Edge*> edge_list;

    Point(double x_, double y_) : x(x_), y(y_), index(-1) {}

    Point(double x_, double y_, int index_) : x(x_), y(y_), index(index_) {}

    void set_zero()
    {
        x = 0.0;
        y = 0.0;
    }

    void set(double x_, double y_)
    {
        x = x_;
        y = y_;
    }

    /// Negate this point.
    Point operator -() const
    {
        Point v;
        v.set(-x, -y);
        return v;
    }

    /// Add a point to this point.
    void operator +=(const Point& v)
    {
        x += v.x;
        y += v.y;
    }

    /// Subtract a point from this point.
    void operator -=(const Point& v)
    {
        x -= v.x;
        y -= v.y;
    }

    /// Multiply this point by a scalar.
    void operator *=(double a)
    {
        x *= a;
        y *= a;
    }

    /// Get the length of this point (the norm).
    double Length() const
    {
        return sqrt(x*x + y*y);
    }

    /// Convert this point into a unit point. Returns the Length.
    double Normalize()
    {
        const double len = Length();
        x /= len;
        y /= len;
        return len;
    }
};

// Represents a simple polygon's edge
struct Edge {

    Point* p, *q;

    /// Constructor
    Edge(Point& p1, Point& p2) : p(&p1), q(&p2)
    {
        if (p1.y > p2.y) {
            q = &p1;
            p = &p2;
        }
        else if (p1.y == p2.y) {
            if (p1.x > p2.x) {
                q = &p1;
                p = &p2;
            }
            else if (p1.x == p2.x) {
                // Repeat points
                assert(false);
            }
        }

        q->edge_list.push_back(this);
    }
};

class Triangle {
public:

    /// Constructor
    Triangle(Point& a, Point& b, Point& c);

    /// Flags to determine if an edge is a Constrained edge
    bool constrained_edge[3];
    /// Flags to determine if an edge is a Delauney edge
    bool delaunay_edge[3];

    Point* GetPoint(int index);
    const Point* GetPoint(int index) const;
    Point* PointCW(const Point& point);
    Point* PointCCW(const Point& point);
    Point* OppositePoint(Triangle& t, const Point& p);

    Triangle* GetNeighbor(int index);
    void MarkNeighbor(Point* p1, Point* p2, Triangle* t);
    void MarkNeighbor(Triangle& t);

    void MarkConstrainedEdge(int index);
    void MarkConstrainedEdge(Edge& edge);
    void MarkConstrainedEdge(Point* p, Point* q);

    int Index(const Point* p);
    int EdgeIndex(const Point* p1, const Point* p2);

    Triangle* NeighborCW(const Point& point);
    Triangle* NeighborCCW(const Point& point);
    bool GetConstrainedEdgeCCW(const Point& p);
    bool GetConstrainedEdgeCW(const Point& p);
    void SetConstrainedEdgeCCW(const Point& p, bool ce);
    void SetConstrainedEdgeCW(const Point& p, bool ce);
    bool GetDelunayEdgeCCW(const Point& p);
    bool GetDelunayEdgeCW(const Point& p);
    void SetDelunayEdgeCCW(const Point& p, bool e);
    void SetDelunayEdgeCW(const Point& p, bool e);

    bool Contains(const Point* p);
    bool Contains(const Edge& e);
    bool Contains(const Point* p, const Point* q);
    void Legalize(Point& point);
    void Legalize(Point& opoint, Point& npoint);
    /**
    * Clears all references to all other triangles and points
    */
    void Clear();
    void ClearNeighbor(const Triangle *triangle);
    void ClearNeighbors();
    void ClearDelunayEdges();

    inline bool IsInterior();
    inline void IsInterior(bool b);

    Triangle& NeighborAcross(const Point& opoint);

    void DebugPrint();

private:

    /// Triangle points
    Point* points_[3];
    /// Neighbor list
    Triangle* neighbors_[3];

    /// Has this triangle been marked as an interior triangle?
    bool interior_;
};

inline bool cmp(const Point* a, const Point* b)
{
    if (a->y < b->y) {
        return true;
    }
    else if (a->y == b->y) {
        // Make sure q is point with greater x value
        if (a->x < b->x) {
            return true;
        }
    }
    return false;
}

/// Add two points_ component-wise.
inline Point operator +(const Point& a, const Point& b)
{
    return Point(a.x + b.x, a.y + b.y);
}

/// Subtract two points_ component-wise.
inline Point operator -(const Point& a, const Point& b)
{
    return Point(a.x - b.x, a.y - b.y);
}

/// Multiply point by scalar
inline Point operator *(double s, const Point& a)
{
    return Point(s * a.x, s * a.y);
}

inline bool operator ==(const Point& a, const Point& b)
{
    return a.x == b.x && a.y == b.y;
}

inline bool operator !=(const Point& a, const Point& b)
{
    return !(a.x == b.x) && !(a.y == b.y);
}

/// Peform the dot product on two vectors.
inline double Dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y;
}

/// Perform the cross product on two vectors. In 2D this produces a scalar.
inline double Cross(const Point& a, const Point& b)
{
    return a.x * b.y - a.y * b.x;
}

/// Perform the cross product on a point and a scalar. In 2D this produces
/// a point.
inline Point Cross(const Point& a, double s)
{
    return Point(s * a.y, -s * a.x);
}

/// Perform the cross product on a scalar and a point. In 2D this produces
/// a point.
inline Point Cross(double s, const Point& a)
{
    return Point(-s * a.y, s * a.x);
}

inline Point* Triangle::GetPoint(int index) 
{
    return points_[index];
}

inline const Point* Triangle::GetPoint(int index) const
{
    return points_[index];
}

inline Triangle* Triangle::GetNeighbor(int index)
{
    return neighbors_[index];
}

inline bool Triangle::Contains(const Point* p)
{
    return p == points_[0] || p == points_[1] || p == points_[2];
}

inline bool Triangle::Contains(const Edge& e)
{
    return Contains(e.p) && Contains(e.q);
}

inline bool Triangle::Contains(const Point* p, const Point* q)
{
    return Contains(p) && Contains(q);
}

inline bool Triangle::IsInterior()
{
    return interior_;
}

inline void Triangle::IsInterior(bool b)
{
    interior_ = b;
}

Triangle::Triangle(Point& a, Point& b, Point& c)
{
    points_[0] = &a; points_[1] = &b; points_[2] = &c;
    neighbors_[0] = NULL; neighbors_[1] = NULL; neighbors_[2] = NULL;
    constrained_edge[0] = constrained_edge[1] = constrained_edge[2] = false;
    delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
    interior_ = false;
}

// Update neighbor pointers
void Triangle::MarkNeighbor(Point* p1, Point* p2, Triangle* t)
{
    if ((p1 == points_[2] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[2]))
        neighbors_[0] = t;
    else if ((p1 == points_[0] && p2 == points_[2]) || (p1 == points_[2] && p2 == points_[0]))
        neighbors_[1] = t;
    else if ((p1 == points_[0] && p2 == points_[1]) || (p1 == points_[1] && p2 == points_[0]))
        neighbors_[2] = t;
    else
        assert(0);
}

// Exhaustive search to update neighbor pointers
void Triangle::MarkNeighbor(Triangle& t)
{
    if (t.Contains(points_[1], points_[2])) {
        neighbors_[0] = &t;
        t.MarkNeighbor(points_[1], points_[2], this);
    }
    else if (t.Contains(points_[0], points_[2])) {
        neighbors_[1] = &t;
        t.MarkNeighbor(points_[0], points_[2], this);
    }
    else if (t.Contains(points_[0], points_[1])) {
        neighbors_[2] = &t;
        t.MarkNeighbor(points_[0], points_[1], this);
    }
}

/**
* Clears all references to all other triangles and points
*/
void Triangle::Clear()
{
    Triangle *t;
    for (int i = 0; i<3; i++)
    {
        t = neighbors_[i];
        if (t != NULL)
        {
            t->ClearNeighbor(this);
        }
    }
    ClearNeighbors();
    points_[0] = points_[1] = points_[2] = NULL;
}

void Triangle::ClearNeighbor(const Triangle *triangle)
{
    if (neighbors_[0] == triangle)
    {
        neighbors_[0] = NULL;
    }
    else if (neighbors_[1] == triangle)
    {
        neighbors_[1] = NULL;
    }
    else
    {
        neighbors_[2] = NULL;
    }
}

void Triangle::ClearNeighbors()
{
    neighbors_[0] = NULL;
    neighbors_[1] = NULL;
    neighbors_[2] = NULL;
}

void Triangle::ClearDelunayEdges()
{
    delaunay_edge[0] = delaunay_edge[1] = delaunay_edge[2] = false;
}

Point* Triangle::OppositePoint(Triangle& t, const Point& p)
{
    Point *cw = t.PointCW(p);
    return PointCW(*cw);
}

// Legalized triangle by rotating clockwise around point(0)
void Triangle::Legalize(Point& point)
{
    points_[1] = points_[0];
    points_[0] = points_[2];
    points_[2] = &point;
}

// Legalize triagnle by rotating clockwise around oPoint
void Triangle::Legalize(Point& opoint, Point& npoint)
{
    if (&opoint == points_[0]) {
        points_[1] = points_[0];
        points_[0] = points_[2];
        points_[2] = &npoint;
    }
    else if (&opoint == points_[1]) {
        points_[2] = points_[1];
        points_[1] = points_[0];
        points_[0] = &npoint;
    }
    else if (&opoint == points_[2]) {
        points_[0] = points_[2];
        points_[2] = points_[1];
        points_[1] = &npoint;
    }
    else {
        assert(0);
    }
}

int Triangle::Index(const Point* p)
{
    if (p == points_[0]) {
        return 0;
    }
    else if (p == points_[1]) {
        return 1;
    }
    else if (p == points_[2]) {
        return 2;
    }
    assert(0);
    return -1;
}

int Triangle::EdgeIndex(const Point* p1, const Point* p2)
{
    if (points_[0] == p1) {
        if (points_[1] == p2) {
            return 2;
        }
        else if (points_[2] == p2) {
            return 1;
        }
    }
    else if (points_[1] == p1) {
        if (points_[2] == p2) {
            return 0;
        }
        else if (points_[0] == p2) {
            return 2;
        }
    }
    else if (points_[2] == p1) {
        if (points_[0] == p2) {
            return 1;
        }
        else if (points_[1] == p2) {
            return 0;
        }
    }
    return -1;
}

void Triangle::MarkConstrainedEdge(int index)
{
    constrained_edge[index] = true;
}

void Triangle::MarkConstrainedEdge(Edge& edge)
{
    MarkConstrainedEdge(edge.p, edge.q);
}

// Mark edge as constrained
void Triangle::MarkConstrainedEdge(Point* p, Point* q)
{
    if ((q == points_[0] && p == points_[1]) || (q == points_[1] && p == points_[0])) {
        constrained_edge[2] = true;
    }
    else if ((q == points_[0] && p == points_[2]) || (q == points_[2] && p == points_[0])) {
        constrained_edge[1] = true;
    }
    else if ((q == points_[1] && p == points_[2]) || (q == points_[2] && p == points_[1])) {
        constrained_edge[0] = true;
    }
}

// The point counter-clockwise to given point
Point* Triangle::PointCW(const Point& point)
{
    if (&point == points_[0]) {
        return points_[2];
    }
    else if (&point == points_[1]) {
        return points_[0];
    }
    else if (&point == points_[2]) {
        return points_[1];
    }
    assert(0);
    return NULL;
}

// The point counter-clockwise to given point
Point* Triangle::PointCCW(const Point& point)
{
    if (&point == points_[0]) {
        return points_[1];
    }
    else if (&point == points_[1]) {
        return points_[2];
    }
    else if (&point == points_[2]) {
        return points_[0];
    }
    assert(0);
    return NULL;
}

// The neighbor clockwise to given point
Triangle* Triangle::NeighborCW(const Point& point)
{
    if (&point == points_[0]) {
        return neighbors_[1];
    }
    else if (&point == points_[1]) {
        return neighbors_[2];
    }
    return neighbors_[0];
}

// The neighbor counter-clockwise to given point
Triangle* Triangle::NeighborCCW(const Point& point)
{
    if (&point == points_[0]) {
        return neighbors_[2];
    }
    else if (&point == points_[1]) {
        return neighbors_[0];
    }
    return neighbors_[1];
}

bool Triangle::GetConstrainedEdgeCCW(const Point& p)
{
    if (&p == points_[0]) {
        return constrained_edge[2];
    }
    else if (&p == points_[1]) {
        return constrained_edge[0];
    }
    return constrained_edge[1];
}

bool Triangle::GetConstrainedEdgeCW(const Point& p)
{
    if (&p == points_[0]) {
        return constrained_edge[1];
    }
    else if (&p == points_[1]) {
        return constrained_edge[2];
    }
    return constrained_edge[0];
}

void Triangle::SetConstrainedEdgeCCW(const Point& p, bool ce)
{
    if (&p == points_[0]) {
        constrained_edge[2] = ce;
    }
    else if (&p == points_[1]) {
        constrained_edge[0] = ce;
    }
    else {
        constrained_edge[1] = ce;
    }
}

void Triangle::SetConstrainedEdgeCW(const Point& p, bool ce)
{
    if (&p == points_[0]) {
        constrained_edge[1] = ce;
    }
    else if (&p == points_[1]) {
        constrained_edge[2] = ce;
    }
    else {
        constrained_edge[0] = ce;
    }
}

bool Triangle::GetDelunayEdgeCCW(const Point& p)
{
    if (&p == points_[0]) {
        return delaunay_edge[2];
    }
    else if (&p == points_[1]) {
        return delaunay_edge[0];
    }
    return delaunay_edge[1];
}

bool Triangle::GetDelunayEdgeCW(const Point& p)
{
    if (&p == points_[0]) {
        return delaunay_edge[1];
    }
    else if (&p == points_[1]) {
        return delaunay_edge[2];
    }
    return delaunay_edge[0];
}

void Triangle::SetDelunayEdgeCCW(const Point& p, bool e)
{
    if (&p == points_[0]) {
        delaunay_edge[2] = e;
    }
    else if (&p == points_[1]) {
        delaunay_edge[0] = e;
    }
    else {
        delaunay_edge[1] = e;
    }
}

void Triangle::SetDelunayEdgeCW(const Point& p, bool e)
{
    if (&p == points_[0]) {
        delaunay_edge[1] = e;
    }
    else if (&p == points_[1]) {
        delaunay_edge[2] = e;
    }
    else {
        delaunay_edge[0] = e;
    }
}

// The neighbor across to given point
Triangle& Triangle::NeighborAcross(const Point& opoint)
{
    if (&opoint == points_[0]) {
        return *neighbors_[0];
    }
    else if (&opoint == points_[1]) {
        return *neighbors_[1];
    }
    return *neighbors_[2];
}

void Triangle::DebugPrint()
{
    using namespace std;
    cout << points_[0]->x << "," << points_[0]->y << " ";
    cout << points_[1]->x << "," << points_[1]->y << " ";
    cout << points_[2]->x << "," << points_[2]->y << endl;
}

const double PI_3div4 = 3 * M_PI / 4;
const double PI_div2 = 1.57079632679489661923;
const double EPSILON = 1e-12;

enum Orientation { CW, CCW, COLLINEAR };

Orientation Orient2d(const Point& pa, const Point& pb, const Point& pc)
{
    double detleft = (pa.x - pc.x) * (pb.y - pc.y);
    double detright = (pa.y - pc.y) * (pb.x - pc.x);
    double val = detleft - detright;
    if (val > -EPSILON && val < EPSILON) {
        return COLLINEAR;
    }
    else if (val > 0) {
        return CCW;
    }
    return CW;
}

bool InScanArea(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
    double oadb = (pa.x - pb.x)*(pd.y - pb.y) - (pd.x - pb.x)*(pa.y - pb.y);
    if (oadb >= -EPSILON) {
        return false;
    }

    double oadc = (pa.x - pc.x)*(pd.y - pc.y) - (pd.x - pc.x)*(pa.y - pc.y);
    if (oadc <= EPSILON) {
        return false;
    }
    return true;
}

struct Node;

// Advancing front node
struct Node {
    Point* point;
    Triangle* triangle;

    Node* next;
    Node* prev;

    double value;

    Node(Point& p) : point(&p), triangle(NULL), next(NULL), prev(NULL), value(p.x)
    {
    }

    Node(Point& p, Triangle& t) : point(&p), triangle(&t), next(NULL), prev(NULL), value(p.x)
    {
    }

};

// Advancing front
class AdvancingFront {
public:

    AdvancingFront(Node& head, Node& tail);
    // Destructor
    ~AdvancingFront();

    Node* head();
    void set_head(Node* node);
    Node* tail();
    void set_tail(Node* node);
    Node* search();
    void set_search(Node* node);

    /// Locate insertion point along advancing front
    Node* LocateNode(double x);

    Node* LocatePoint(const Point* point);

private:

    Node* head_, *tail_, *search_node_;

    Node* FindSearchNode(double x);
};

inline Node* AdvancingFront::head()
{
    return head_;
}
inline void AdvancingFront::set_head(Node* node)
{
    head_ = node;
}

inline Node* AdvancingFront::tail()
{
    return tail_;
}
inline void AdvancingFront::set_tail(Node* node)
{
    tail_ = node;
}

inline Node* AdvancingFront::search()
{
    return search_node_;
}

inline void AdvancingFront::set_search(Node* node)
{
    search_node_ = node;
}

AdvancingFront::AdvancingFront(Node& head, Node& tail)
{
    head_ = &head;
    tail_ = &tail;
    search_node_ = &head;
}

Node* AdvancingFront::LocateNode(double x)
{
    Node* node = search_node_;

    if (x < node->value) {
        while ((node = node->prev) != NULL) {
            if (x >= node->value) {
                search_node_ = node;
                return node;
            }
        }
    }
    else {
        while ((node = node->next) != NULL) {
            if (x < node->value) {
                search_node_ = node->prev;
                return node->prev;
            }
        }
    }
    return NULL;
}

Node* AdvancingFront::FindSearchNode(double x)
{
    (void)x; // suppress compiler warnings "unused parameter 'x'"
    // TODO: implement BST index
    return search_node_;
}

Node* AdvancingFront::LocatePoint(const Point* point)
{
    const double px = point->x;
    Node* node = FindSearchNode(px);
    const double nx = node->point->x;

    if (px == nx) {
        if (point != node->point) {
            // We might have two nodes with same x value for a short time
            if (point == node->prev->point) {
                node = node->prev;
            }
            else if (point == node->next->point) {
                node = node->next;
            }
            else {
                assert(0);
            }
        }
    }
    else if (px < nx) {
        while ((node = node->prev) != NULL) {
            if (point == node->point) {
                break;
            }
        }
    }
    else {
        while ((node = node->next) != NULL) {
            if (point == node->point)
                break;
        }
    }
    if (node) search_node_ = node;
    return node;
}

AdvancingFront::~AdvancingFront()
{
}

const double kAlpha = 0.3;

struct Point;
class Triangle;
struct Node;
struct Edge;
class AdvancingFront;

class SweepContext {
public:

    /// Constructor
    SweepContext(const std::vector<Point*>& polyline);
    /// Destructor
    ~SweepContext();

    void set_head(Point* p1);

    Point* head() const;

    void set_tail(Point* p1);

    Point* tail() const;

    size_t point_count() const;

    Node& LocateNode(const Point& point);

    void RemoveNode(Node* node);

    void CreateAdvancingFront(const std::vector<Node*>& nodes);

    /// Try to map a node to all sides of this triangle that don't have a neighbor
    void MapTriangleToNodes(Triangle& t);

    void AddToMap(Triangle* triangle);

    Point* GetPoint(size_t index);

    Point* GetPoints();

    void RemoveFromMap(Triangle* triangle);

    void AddHole(const std::vector<Point*>& polyline);

    void AddPoint(Point* point);

    AdvancingFront* front() const;

    void MeshClean(Triangle& triangle);

    std::vector<Triangle*> &GetTriangles();
    std::list<Triangle*> &GetMap();

    std::vector<Edge*> edge_list;

    struct Basin {
        Node* left_node;
        Node* bottom_node;
        Node* right_node;
        double width;
        bool left_highest;

        Basin() : left_node(NULL), bottom_node(NULL), right_node(NULL), width(0.0), left_highest(false)
        {
        }

        void Clear()
        {
            left_node = NULL;
            bottom_node = NULL;
            right_node = NULL;
            width = 0.0;
            left_highest = false;
        }
    };

    struct EdgeEvent {
        Edge* constrained_edge;
        bool right;

        EdgeEvent() : constrained_edge(NULL), right(false)
        {
        }
    };

    Basin basin;
    EdgeEvent edge_event;

private:

    friend class Sweep;

    std::vector<Triangle*> triangles_;
    std::list<Triangle*> map_;
    std::vector<Point*> points_;

    // Advancing front
    AdvancingFront* front_;
    // head point used with advancing front
    Point* head_;
    // tail point used with advancing front
    Point* tail_;

    Node *af_head_, *af_middle_, *af_tail_;

    void InitTriangulation();
    void InitEdges(const std::vector<Point*>& polyline);

};

inline AdvancingFront* SweepContext::front() const
{
    return front_;
}

inline size_t SweepContext::point_count() const
{
    return points_.size();
}

inline void SweepContext::set_head(Point* p1)
{
    head_ = p1;
}

inline Point* SweepContext::head() const
{
    return head_;
}

inline void SweepContext::set_tail(Point* p1)
{
    tail_ = p1;
}

inline Point* SweepContext::tail() const
{
    return tail_;
}

SweepContext::SweepContext(const std::vector<Point*>& polyline) 
: points_(polyline)
, front_(0)
, head_(0)
, tail_(0)
, af_head_(0)
, af_middle_(0)
, af_tail_(0)
{
    InitEdges(points_);
}

void SweepContext::AddHole(const std::vector<Point*>& polyline)
{
    InitEdges(polyline);
    for (unsigned int i = 0; i < polyline.size(); i++) {
        points_.push_back(polyline[i]);
    }
}

void SweepContext::AddPoint(Point* point) {
    points_.push_back(point);
}

std::vector<Triangle*> &SweepContext::GetTriangles()
{
    return triangles_;
}

std::list<Triangle*> &SweepContext::GetMap()
{
    return map_;
}

void SweepContext::InitTriangulation()
{
    double xmax(points_[0]->x), xmin(points_[0]->x);
    double ymax(points_[0]->y), ymin(points_[0]->y);

    // Calculate bounds.
    for (unsigned int i = 0; i < points_.size(); i++) {
        Point& p = *points_[i];
        if (p.x > xmax)
            xmax = p.x;
        if (p.x < xmin)
            xmin = p.x;
        if (p.y > ymax)
            ymax = p.y;
        if (p.y < ymin)
            ymin = p.y;
    }

    double dx = kAlpha * (xmax - xmin);
    double dy = kAlpha * (ymax - ymin);
    head_ = new Point(xmax + dx, ymin - dy);
    tail_ = new Point(xmin - dx, ymin - dy);

    // Sort points along y-axis
    std::sort(points_.begin(), points_.end(), cmp);

}

void SweepContext::InitEdges(const std::vector<Point*>& polyline)
{
    size_t num_points = polyline.size();
    if (num_points > 1) {
        for (size_t i = 0; i < num_points; i++) {
            size_t j = i < num_points - 1 ? i + 1 : 0;
            edge_list.push_back(new Edge(*polyline[i], *polyline[j]));
        }
    }
}

Point* SweepContext::GetPoint(size_t index)
{
    return points_[index];
}

void SweepContext::AddToMap(Triangle* triangle)
{
    map_.push_back(triangle);
}

Node& SweepContext::LocateNode(const Point& point)
{
    // TODO implement search tree
    return *front_->LocateNode(point.x);
}

void SweepContext::CreateAdvancingFront(const std::vector<Node*>& nodes)
{

    (void)nodes;
    // Initial triangle
    Triangle* triangle = new Triangle(*points_[0], *tail_, *head_);

    map_.push_back(triangle);

    af_head_ = new Node(*triangle->GetPoint(1), *triangle);
    af_middle_ = new Node(*triangle->GetPoint(0), *triangle);
    af_tail_ = new Node(*triangle->GetPoint(2));
    front_ = new AdvancingFront(*af_head_, *af_tail_);

    // TODO: More intuitive if head is middles next and not previous?
    //       so swap head and tail
    af_head_->next = af_middle_;
    af_middle_->next = af_tail_;
    af_middle_->prev = af_head_;
    af_tail_->prev = af_middle_;
}

void SweepContext::RemoveNode(Node* node)
{
    delete node;
}

void SweepContext::MapTriangleToNodes(Triangle& t)
{
    for (int i = 0; i < 3; i++) {
        if (!t.GetNeighbor(i)) {
            Node* n = front_->LocatePoint(t.PointCW(*t.GetPoint(i)));
            if (n)
                n->triangle = &t;
        }
    }
}

void SweepContext::RemoveFromMap(Triangle* triangle)
{
    map_.remove(triangle);
}

void SweepContext::MeshClean(Triangle& triangle)
{
    std::vector<Triangle *> triangles;
    triangles.push_back(&triangle);

    while (!triangles.empty()){
        Triangle *t = triangles.back();
        triangles.pop_back();

        if (t != NULL && !t->IsInterior()) {
            t->IsInterior(true);
            triangles_.push_back(t);
            for (int i = 0; i < 3; i++) {
                if (!t->constrained_edge[i])
                    triangles.push_back(t->GetNeighbor(i));
            }
        }
    }
}

SweepContext::~SweepContext()
{

    // Clean up memory

    delete head_;
    delete tail_;
    delete front_;
    delete af_head_;
    delete af_middle_;
    delete af_tail_;

    typedef std::list<Triangle*> type_list;

    for (type_list::iterator iter = map_.begin(); iter != map_.end(); ++iter) {
        Triangle* ptr = *iter;
        delete ptr;
    }

    for (unsigned int i = 0; i < edge_list.size(); i++) {
        delete edge_list[i];
    }

}

class Sweep
{
public:

    void Triangulate(SweepContext& tcx);

    ~Sweep();

private:

    void SweepPoints(SweepContext& tcx);

    Node& PointEvent(SweepContext& tcx, Point& point);

    void EdgeEvent(SweepContext& tcx, Edge* edge, Node* node);

    void EdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle* triangle, Point& point);

    Node& NewFrontTriangle(SweepContext& tcx, Point& point, Node& node);

    void Fill(SweepContext& tcx, Node& node);

    bool Legalize(SweepContext& tcx, Triangle& t);

    bool Incircle(const Point& pa, const Point& pb, const Point& pc, const Point& pd) const;

    /**
    * Rotates a triangle pair one vertex CW
    *<pre>
    *       n2                    n2
    *  P +-----+             P +-----+
    *    | t  /|               |\  t |
    *    |   / |               | \   |
    *  n1|  /  |n3           n1|  \  |n3
    *    | /   |    after CW   |   \ |
    *    |/ oT |               | oT \|
    *    +-----+ oP            +-----+
    *       n4                    n4
    * </pre>
    */
    void RotateTrianglePair(Triangle& t, Point& p, Triangle& ot, Point& op) const;

    /**
    * Fills holes in the Advancing Front
    *
    *
    * @param tcx
    * @param n
    */
    void FillAdvancingFront(SweepContext& tcx, Node& n);

    // Decision-making about when to Fill hole.
    // Contributed by ToolmakerSteve2
    bool LargeHole_DontFill(const Node* node) const;
    bool AngleExceeds90Degrees(const Point* origin, const Point* pa, const Point* pb) const;
    bool AngleExceedsPlus90DegreesOrIsNegative(const Point* origin, const Point* pa, const Point* pb) const;
    double Angle(const Point* origin, const Point* pa, const Point* pb) const;

    /**
    *
    * @param node - middle node
    * @return the angle between 3 front nodes
    */
    double HoleAngle(const Node& node) const;

    /**
    * The basin angle is decided against the horizontal line [1,0]
    */
    double BasinAngle(const Node& node) const;

    /**
    * Fills a basin that has formed on the Advancing Front to the right
    * of given node.<br>
    * First we decide a left,bottom and right node that forms the
    * boundaries of the basin. Then we do a reqursive fill.
    *
    * @param tcx
    * @param node - starting node, this or next node will be left node
    */
    void FillBasin(SweepContext& tcx, Node& node);

    /**
    * Recursive algorithm to fill a Basin with triangles
    *
    * @param tcx
    * @param node - bottom_node
    * @param cnt - counter used to alternate on even and odd numbers
    */
    void FillBasinReq(SweepContext& tcx, Node* node);

    bool IsShallow(SweepContext& tcx, Node& node);

    bool IsEdgeSideOfTriangle(Triangle& triangle, Point& ep, Point& eq);

    void FillEdgeEvent(SweepContext& tcx, Edge* edge, Node* node);

    void FillRightAboveEdgeEvent(SweepContext& tcx, Edge* edge, Node* node);

    void FillRightBelowEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FillRightConcaveEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FillRightConvexEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FillLeftAboveEdgeEvent(SweepContext& tcx, Edge* edge, Node* node);

    void FillLeftBelowEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FillLeftConcaveEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FillLeftConvexEdgeEvent(SweepContext& tcx, Edge* edge, Node& node);

    void FlipEdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle* t, Point& p);

    /**
    * After a flip we have two triangles and know that only one will still be
    * intersecting the edge. So decide which to contiune with and legalize the other
    *
    * @param tcx
    * @param o - should be the result of an orient2d( eq, op, ep )
    * @param t - triangle 1
    * @param ot - triangle 2
    * @param p - a point shared by both triangles
    * @param op - another point shared by both triangles
    * @return returns the triangle still intersecting the edge
    */
    Triangle& NextFlipTriangle(SweepContext& tcx, int o, Triangle&  t, Triangle& ot, Point& p, Point& op);

    /**
    * When we need to traverse from one triangle to the next we need
    * the point in current triangle that is the opposite point to the next
    * triangle.
    *
    * @param ep
    * @param eq
    * @param ot
    * @param op
    * @return
    */
    Point& NextFlipPoint(Point& ep, Point& eq, Triangle& ot, Point& op);

    /**
    * Scan part of the FlipScan algorithm<br>
    * When a triangle pair isn't flippable we will scan for the next
    * point that is inside the flip triangle scan area. When found
    * we generate a new flipEdgeEvent
    *
    * @param tcx
    * @param ep - last point on the edge we are traversing
    * @param eq - first point on the edge we are traversing
    * @param flipTriangle - the current triangle sharing the point eq with edge
    * @param t
    * @param p
    */
    void FlipScanEdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle& flip_triangle, Triangle& t, Point& p);

    void FinalizationPolygon(SweepContext& tcx);

    std::vector<Node*> nodes_;

};


// Triangulate simple polygon with holes
void Sweep::Triangulate(SweepContext& tcx)
{
    tcx.InitTriangulation();
    tcx.CreateAdvancingFront(nodes_);
    // Sweep points; build mesh
    SweepPoints(tcx);
    // Clean up
    FinalizationPolygon(tcx);
}

void Sweep::SweepPoints(SweepContext& tcx)
{
    for (size_t i = 1; i < tcx.point_count(); i++) {
        Point& point = *tcx.GetPoint(i);
        Node* node = &PointEvent(tcx, point);
        for (unsigned int i = 0; i < point.edge_list.size(); i++) {
            EdgeEvent(tcx, point.edge_list[i], node);
        }
    }
}

void Sweep::FinalizationPolygon(SweepContext& tcx)
{
    // Get an Internal triangle to start with
    Triangle* t = tcx.front()->head()->next->triangle;
    Point* p = tcx.front()->head()->next->point;
    while (!t->GetConstrainedEdgeCW(*p)) {
        t = t->NeighborCCW(*p);
    }

    // Collect interior triangles constrained by edges
    tcx.MeshClean(*t);
}

Node& Sweep::PointEvent(SweepContext& tcx, Point& point)
{
    Node& node = tcx.LocateNode(point);
    Node& new_node = NewFrontTriangle(tcx, point, node);

    // Only need to check +epsilon since point never have smaller
    // x value than node due to how we fetch nodes from the front
    if (point.x <= node.point->x + EPSILON) {
        Fill(tcx, node);
    }

    //tcx.AddNode(new_node);

    FillAdvancingFront(tcx, new_node);
    return new_node;
}

void Sweep::EdgeEvent(SweepContext& tcx, Edge* edge, Node* node)
{
    tcx.edge_event.constrained_edge = edge;
    tcx.edge_event.right = (edge->p->x > edge->q->x);

    if (IsEdgeSideOfTriangle(*node->triangle, *edge->p, *edge->q)) {
        return;
    }

    // For now we will do all needed filling
    // TODO: integrate with flip process might give some better performance
    //       but for now this avoid the issue with cases that needs both flips and fills
    FillEdgeEvent(tcx, edge, node);
    EdgeEvent(tcx, *edge->p, *edge->q, node->triangle, *edge->q);
}

void Sweep::EdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle* triangle, Point& point)
{
    if (IsEdgeSideOfTriangle(*triangle, ep, eq)) {
        return;
    }

    Point* p1 = triangle->PointCCW(point);
    Orientation o1 = Orient2d(eq, *p1, ep);
    if (o1 == COLLINEAR) {
        if (triangle->Contains(&eq, p1)) {
            triangle->MarkConstrainedEdge(&eq, p1);
            // We are modifying the constraint maybe it would be better to
            // not change the given constraint and just keep a variable for the new constraint
            tcx.edge_event.constrained_edge->q = p1;
            triangle = &triangle->NeighborAcross(point);
            EdgeEvent(tcx, ep, *p1, triangle, *p1);
        }
        else {
            std::runtime_error("EdgeEvent - collinear points not supported");
            assert(0);
        }
        return;
    }

    Point* p2 = triangle->PointCW(point);
    Orientation o2 = Orient2d(eq, *p2, ep);
    if (o2 == COLLINEAR) {
        if (triangle->Contains(&eq, p2)) {
            triangle->MarkConstrainedEdge(&eq, p2);
            // We are modifying the constraint maybe it would be better to
            // not change the given constraint and just keep a variable for the new constraint
            tcx.edge_event.constrained_edge->q = p2;
            triangle = &triangle->NeighborAcross(point);
            EdgeEvent(tcx, ep, *p2, triangle, *p2);
        }
        else {
            std::runtime_error("EdgeEvent - collinear points not supported");
            assert(0);
        }
        return;
    }

    if (o1 == o2) {
        // Need to decide if we are rotating CW or CCW to get to a triangle
        // that will cross edge
        if (o1 == CW) {
            triangle = triangle->NeighborCCW(point);
        }
        else{
            triangle = triangle->NeighborCW(point);
        }
        EdgeEvent(tcx, ep, eq, triangle, point);
    }
    else {
        // This triangle crosses constraint so lets flippin start!
        FlipEdgeEvent(tcx, ep, eq, triangle, point);
    }
}

bool Sweep::IsEdgeSideOfTriangle(Triangle& triangle, Point& ep, Point& eq)
{
    const int index = triangle.EdgeIndex(&ep, &eq);

    if (index != -1) {
        triangle.MarkConstrainedEdge(index);
        Triangle* t = triangle.GetNeighbor(index);
        if (t) {
            t->MarkConstrainedEdge(&ep, &eq);
        }
        return true;
    }
    return false;
}

Node& Sweep::NewFrontTriangle(SweepContext& tcx, Point& point, Node& node)
{
    Triangle* triangle = new Triangle(point, *node.point, *node.next->point);

    triangle->MarkNeighbor(*node.triangle);
    tcx.AddToMap(triangle);

    Node* new_node = new Node(point);
    nodes_.push_back(new_node);

    new_node->next = node.next;
    new_node->prev = &node;
    node.next->prev = new_node;
    node.next = new_node;

    if (!Legalize(tcx, *triangle)) {
        tcx.MapTriangleToNodes(*triangle);
    }

    return *new_node;
}

void Sweep::Fill(SweepContext& tcx, Node& node)
{
    Triangle* triangle = new Triangle(*node.prev->point, *node.point, *node.next->point);

    // TODO: should copy the constrained_edge value from neighbor triangles
    //       for now constrained_edge values are copied during the legalize
    triangle->MarkNeighbor(*node.prev->triangle);
    triangle->MarkNeighbor(*node.triangle);

    tcx.AddToMap(triangle);

    // Update the advancing front
    node.prev->next = node.next;
    node.next->prev = node.prev;

    // If it was legalized the triangle has already been mapped
    if (!Legalize(tcx, *triangle)) {
        tcx.MapTriangleToNodes(*triangle);
    }

}

void Sweep::FillAdvancingFront(SweepContext& tcx, Node& n)
{

    // Fill right holes
    Node* node = n.next;

    while (node->next) {
        // if HoleAngle exceeds 90 degrees then break.
        if (LargeHole_DontFill(node)) break;
        Fill(tcx, *node);
        node = node->next;
    }

    // Fill left holes
    node = n.prev;

    while (node->prev) {
        // if HoleAngle exceeds 90 degrees then break.
        if (LargeHole_DontFill(node)) break;
        Fill(tcx, *node);
        node = node->prev;
    }

    // Fill right basins
    if (n.next && n.next->next) {
        const double angle = BasinAngle(n);
        if (angle < PI_3div4) {
            FillBasin(tcx, n);
        }
    }
}

// True if HoleAngle exceeds 90 degrees.
bool Sweep::LargeHole_DontFill(const Node* node) const {

    const Node* nextNode = node->next;
    const Node* prevNode = node->prev;
    if (!AngleExceeds90Degrees(node->point, nextNode->point, prevNode->point))
        return false;

    // Check additional points on front.
    const Node* next2Node = nextNode->next;
    // "..Plus.." because only want angles on same side as point being added.
    if ((next2Node != NULL) && !AngleExceedsPlus90DegreesOrIsNegative(node->point, next2Node->point, prevNode->point))
        return false;

    const Node* prev2Node = prevNode->prev;
    // "..Plus.." because only want angles on same side as point being added.
    if ((prev2Node != NULL) && !AngleExceedsPlus90DegreesOrIsNegative(node->point, nextNode->point, prev2Node->point))
        return false;

    return true;
}

bool Sweep::AngleExceeds90Degrees(const Point* origin, const Point* pa, const Point* pb) const {
    const double angle = Angle(origin, pa, pb);
    return ((angle > PI_div2) || (angle < -PI_div2));
}

bool Sweep::AngleExceedsPlus90DegreesOrIsNegative(const Point* origin, const Point* pa, const Point* pb) const {
    const double angle = Angle(origin, pa, pb);
    return (angle > PI_div2) || (angle < 0);
}

double Sweep::Angle(const Point* origin, const Point* pa, const Point* pb) const {
    /* Complex plane
    * ab = cosA +i*sinA
    * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
    * atan2(y,x) computes the principal value of the argument function
    * applied to the complex number x+iy
    * Where x = ax*bx + ay*by
    *       y = ax*by - ay*bx
    */
    const double px = origin->x;
    const double py = origin->y;
    const double ax = pa->x - px;
    const double ay = pa->y - py;
    const double bx = pb->x - px;
    const double by = pb->y - py;
    const double x = ax * by - ay * bx;
    const double y = ax * bx + ay * by;
    return atan2(x, y);
}

double Sweep::BasinAngle(const Node& node) const
{
    const double ax = node.point->x - node.next->next->point->x;
    const double ay = node.point->y - node.next->next->point->y;
    return atan2(ay, ax);
}

double Sweep::HoleAngle(const Node& node) const
{
    /* Complex plane
    * ab = cosA +i*sinA
    * ab = (ax + ay*i)(bx + by*i) = (ax*bx + ay*by) + i(ax*by-ay*bx)
    * atan2(y,x) computes the principal value of the argument function
    * applied to the complex number x+iy
    * Where x = ax*bx + ay*by
    *       y = ax*by - ay*bx
    */
    const double ax = node.next->point->x - node.point->x;
    const double ay = node.next->point->y - node.point->y;
    const double bx = node.prev->point->x - node.point->x;
    const double by = node.prev->point->y - node.point->y;
    return atan2(ax * by - ay * bx, ax * bx + ay * by);
}

bool Sweep::Legalize(SweepContext& tcx, Triangle& t)
{
    // To legalize a triangle we start by finding if any of the three edges
    // violate the Delaunay condition
    for (int i = 0; i < 3; i++) {
        if (t.delaunay_edge[i])
            continue;

        Triangle* ot = t.GetNeighbor(i);

        if (ot) {
            Point* p = t.GetPoint(i);
            Point* op = ot->OppositePoint(t, *p);
            int oi = ot->Index(op);

            // If this is a Constrained Edge or a Delaunay Edge(only during recursive legalization)
            // then we should not try to legalize
            if (ot->constrained_edge[oi] || ot->delaunay_edge[oi]) {
                t.constrained_edge[i] = ot->constrained_edge[oi];
                continue;
            }

            bool inside = Incircle(*p, *t.PointCCW(*p), *t.PointCW(*p), *op);

            if (inside) {
                // Lets mark this shared edge as Delaunay
                t.delaunay_edge[i] = true;
                ot->delaunay_edge[oi] = true;

                // Lets rotate shared edge one vertex CW to legalize it
                RotateTrianglePair(t, *p, *ot, *op);

                // We now got one valid Delaunay Edge shared by two triangles
                // This gives us 4 new edges to check for Delaunay

                // Make sure that triangle to node mapping is done only one time for a specific triangle
                bool not_legalized = !Legalize(tcx, t);
                if (not_legalized) {
                    tcx.MapTriangleToNodes(t);
                }

                not_legalized = !Legalize(tcx, *ot);
                if (not_legalized)
                    tcx.MapTriangleToNodes(*ot);

                // Reset the Delaunay edges, since they only are valid Delaunay edges
                // until we add a new triangle or point.
                // XXX: need to think about this. Can these edges be tried after we
                //      return to previous recursive level?
                t.delaunay_edge[i] = false;
                ot->delaunay_edge[oi] = false;

                // If triangle have been legalized no need to check the other edges since
                // the recursive legalization will handles those so we can end here.
                return true;
            }
        }
    }
    return false;
}

bool Sweep::Incircle(const Point& pa, const Point& pb, const Point& pc, const Point& pd) const
{
    const double adx = pa.x - pd.x;
    const double ady = pa.y - pd.y;
    const double bdx = pb.x - pd.x;
    const double bdy = pb.y - pd.y;

    const double adxbdy = adx * bdy;
    const double bdxady = bdx * ady;
    const double oabd = adxbdy - bdxady;

    if (oabd <= 0)
        return false;

    const double cdx = pc.x - pd.x;
    const double cdy = pc.y - pd.y;

    const double cdxady = cdx * ady;
    const double adxcdy = adx * cdy;
    const double ocad = cdxady - adxcdy;

    if (ocad <= 0)
        return false;

    const double bdxcdy = bdx * cdy;
    const double cdxbdy = cdx * bdy;

    const double alift = adx * adx + ady * ady;
    const double blift = bdx * bdx + bdy * bdy;
    const double clift = cdx * cdx + cdy * cdy;

    const double det = alift * (bdxcdy - cdxbdy) + blift * ocad + clift * oabd;

    return det > 0;
}

void Sweep::RotateTrianglePair(Triangle& t, Point& p, Triangle& ot, Point& op) const
{
    Triangle* n1, *n2, *n3, *n4;
    n1 = t.NeighborCCW(p);
    n2 = t.NeighborCW(p);
    n3 = ot.NeighborCCW(op);
    n4 = ot.NeighborCW(op);

    bool ce1, ce2, ce3, ce4;
    ce1 = t.GetConstrainedEdgeCCW(p);
    ce2 = t.GetConstrainedEdgeCW(p);
    ce3 = ot.GetConstrainedEdgeCCW(op);
    ce4 = ot.GetConstrainedEdgeCW(op);

    bool de1, de2, de3, de4;
    de1 = t.GetDelunayEdgeCCW(p);
    de2 = t.GetDelunayEdgeCW(p);
    de3 = ot.GetDelunayEdgeCCW(op);
    de4 = ot.GetDelunayEdgeCW(op);

    t.Legalize(p, op);
    ot.Legalize(op, p);

    // Remap delaunay_edge
    ot.SetDelunayEdgeCCW(p, de1);
    t.SetDelunayEdgeCW(p, de2);
    t.SetDelunayEdgeCCW(op, de3);
    ot.SetDelunayEdgeCW(op, de4);

    // Remap constrained_edge
    ot.SetConstrainedEdgeCCW(p, ce1);
    t.SetConstrainedEdgeCW(p, ce2);
    t.SetConstrainedEdgeCCW(op, ce3);
    ot.SetConstrainedEdgeCW(op, ce4);

    // Remap neighbors
    // XXX: might optimize the markNeighbor by keeping track of
    //      what side should be assigned to what neighbor after the
    //      rotation. Now mark neighbor does lots of testing to find
    //      the right side.
    t.ClearNeighbors();
    ot.ClearNeighbors();
    if (n1) ot.MarkNeighbor(*n1);
    if (n2) t.MarkNeighbor(*n2);
    if (n3) t.MarkNeighbor(*n3);
    if (n4) ot.MarkNeighbor(*n4);
    t.MarkNeighbor(ot);
}

void Sweep::FillBasin(SweepContext& tcx, Node& node)
{
    if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
        tcx.basin.left_node = node.next->next;
    }
    else {
        tcx.basin.left_node = node.next;
    }

    // Find the bottom and right node
    tcx.basin.bottom_node = tcx.basin.left_node;
    while (tcx.basin.bottom_node->next
        && tcx.basin.bottom_node->point->y >= tcx.basin.bottom_node->next->point->y) {
        tcx.basin.bottom_node = tcx.basin.bottom_node->next;
    }
    if (tcx.basin.bottom_node == tcx.basin.left_node) {
        // No valid basin
        return;
    }

    tcx.basin.right_node = tcx.basin.bottom_node;
    while (tcx.basin.right_node->next
        && tcx.basin.right_node->point->y < tcx.basin.right_node->next->point->y) {
        tcx.basin.right_node = tcx.basin.right_node->next;
    }
    if (tcx.basin.right_node == tcx.basin.bottom_node) {
        // No valid basins
        return;
    }

    tcx.basin.width = tcx.basin.right_node->point->x - tcx.basin.left_node->point->x;
    tcx.basin.left_highest = tcx.basin.left_node->point->y > tcx.basin.right_node->point->y;

    FillBasinReq(tcx, tcx.basin.bottom_node);
}

void Sweep::FillBasinReq(SweepContext& tcx, Node* node)
{
    // if shallow stop filling
    if (IsShallow(tcx, *node)) {
        return;
    }

    Fill(tcx, *node);

    if (node->prev == tcx.basin.left_node && node->next == tcx.basin.right_node) {
        return;
    }
    else if (node->prev == tcx.basin.left_node) {
        Orientation o = Orient2d(*node->point, *node->next->point, *node->next->next->point);
        if (o == CW) {
            return;
        }
        node = node->next;
    }
    else if (node->next == tcx.basin.right_node) {
        Orientation o = Orient2d(*node->point, *node->prev->point, *node->prev->prev->point);
        if (o == CCW) {
            return;
        }
        node = node->prev;
    }
    else {
        // Continue with the neighbor node with lowest Y value
        if (node->prev->point->y < node->next->point->y) {
            node = node->prev;
        }
        else {
            node = node->next;
        }
    }

    FillBasinReq(tcx, node);
}

bool Sweep::IsShallow(SweepContext& tcx, Node& node)
{
    double height;

    if (tcx.basin.left_highest) {
        height = tcx.basin.left_node->point->y - node.point->y;
    }
    else {
        height = tcx.basin.right_node->point->y - node.point->y;
    }

    // if shallow stop filling
    if (tcx.basin.width > height) {
        return true;
    }
    return false;
}

void Sweep::FillEdgeEvent(SweepContext& tcx, Edge* edge, Node* node)
{
    if (tcx.edge_event.right) {
        FillRightAboveEdgeEvent(tcx, edge, node);
    }
    else {
        FillLeftAboveEdgeEvent(tcx, edge, node);
    }
}

void Sweep::FillRightAboveEdgeEvent(SweepContext& tcx, Edge* edge, Node* node)
{
    while (node->next->point->x < edge->p->x) {
        // Check if next node is below the edge
        if (Orient2d(*edge->q, *node->next->point, *edge->p) == CCW) {
            FillRightBelowEdgeEvent(tcx, edge, *node);
        }
        else {
            node = node->next;
        }
    }
}

void Sweep::FillRightBelowEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    if (node.point->x < edge->p->x) {
        if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
            // Concave
            FillRightConcaveEdgeEvent(tcx, edge, node);
        }
        else{
            // Convex
            FillRightConvexEdgeEvent(tcx, edge, node);
            // Retry this one
            FillRightBelowEdgeEvent(tcx, edge, node);
        }
    }
}

void Sweep::FillRightConcaveEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    Fill(tcx, *node.next);
    if (node.next->point != edge->p) {
        // Next above or below edge?
        if (Orient2d(*edge->q, *node.next->point, *edge->p) == CCW) {
            // Below
            if (Orient2d(*node.point, *node.next->point, *node.next->next->point) == CCW) {
                // Next is concave
                FillRightConcaveEdgeEvent(tcx, edge, node);
            }
            else {
                // Next is convex
            }
        }
    }

}

void Sweep::FillRightConvexEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    // Next concave or convex?
    if (Orient2d(*node.next->point, *node.next->next->point, *node.next->next->next->point) == CCW) {
        // Concave
        FillRightConcaveEdgeEvent(tcx, edge, *node.next);
    }
    else{
        // Convex
        // Next above or below edge?
        if (Orient2d(*edge->q, *node.next->next->point, *edge->p) == CCW) {
            // Below
            FillRightConvexEdgeEvent(tcx, edge, *node.next);
        }
        else{
            // Above
        }
    }
}

void Sweep::FillLeftAboveEdgeEvent(SweepContext& tcx, Edge* edge, Node* node)
{
    while (node->prev->point->x > edge->p->x) {
        // Check if next node is below the edge
        if (Orient2d(*edge->q, *node->prev->point, *edge->p) == CW) {
            FillLeftBelowEdgeEvent(tcx, edge, *node);
        }
        else {
            node = node->prev;
        }
    }
}

void Sweep::FillLeftBelowEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    if (node.point->x > edge->p->x) {
        if (Orient2d(*node.point, *node.prev->point, *node.prev->prev->point) == CW) {
            // Concave
            FillLeftConcaveEdgeEvent(tcx, edge, node);
        }
        else {
            // Convex
            FillLeftConvexEdgeEvent(tcx, edge, node);
            // Retry this one
            FillLeftBelowEdgeEvent(tcx, edge, node);
        }
    }
}

void Sweep::FillLeftConvexEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    // Next concave or convex?
    if (Orient2d(*node.prev->point, *node.prev->prev->point, *node.prev->prev->prev->point) == CW) {
        // Concave
        FillLeftConcaveEdgeEvent(tcx, edge, *node.prev);
    }
    else{
        // Convex
        // Next above or below edge?
        if (Orient2d(*edge->q, *node.prev->prev->point, *edge->p) == CW) {
            // Below
            FillLeftConvexEdgeEvent(tcx, edge, *node.prev);
        }
        else{
            // Above
        }
    }
}

void Sweep::FillLeftConcaveEdgeEvent(SweepContext& tcx, Edge* edge, Node& node)
{
    Fill(tcx, *node.prev);
    if (node.prev->point != edge->p) {
        // Next above or below edge?
        if (Orient2d(*edge->q, *node.prev->point, *edge->p) == CW) {
            // Below
            if (Orient2d(*node.point, *node.prev->point, *node.prev->prev->point) == CW) {
                // Next is concave
                FillLeftConcaveEdgeEvent(tcx, edge, node);
            }
            else{
                // Next is convex
            }
        }
    }

}

void Sweep::FlipEdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle* t, Point& p)
{
    Triangle& ot = t->NeighborAcross(p);
    Point& op = *ot.OppositePoint(*t, p);

    if (InScanArea(p, *t->PointCCW(p), *t->PointCW(p), op)) {
        // Lets rotate shared edge one vertex CW
        RotateTrianglePair(*t, p, ot, op);
        tcx.MapTriangleToNodes(*t);
        tcx.MapTriangleToNodes(ot);

        if (p == eq && op == ep) {
            if (eq == *tcx.edge_event.constrained_edge->q && ep == *tcx.edge_event.constrained_edge->p) {
                t->MarkConstrainedEdge(&ep, &eq);
                ot.MarkConstrainedEdge(&ep, &eq);
                Legalize(tcx, *t);
                Legalize(tcx, ot);
            }
            else {
                // XXX: I think one of the triangles should be legalized here?
            }
        }
        else {
            Orientation o = Orient2d(eq, op, ep);
            t = &NextFlipTriangle(tcx, (int)o, *t, ot, p, op);
            FlipEdgeEvent(tcx, ep, eq, t, p);
        }
    }
    else {
        Point& newP = NextFlipPoint(ep, eq, ot, op);
        FlipScanEdgeEvent(tcx, ep, eq, *t, ot, newP);
        EdgeEvent(tcx, ep, eq, t, p);
    }
}

Triangle& Sweep::NextFlipTriangle(SweepContext& tcx, int o, Triangle& t, Triangle& ot, Point& p, Point& op)
{
    if (o == CCW) {
        // ot is not crossing edge after flip
        int edge_index = ot.EdgeIndex(&p, &op);
        ot.delaunay_edge[edge_index] = true;
        Legalize(tcx, ot);
        ot.ClearDelunayEdges();
        return t;
    }

    // t is not crossing edge after flip
    int edge_index = t.EdgeIndex(&p, &op);

    t.delaunay_edge[edge_index] = true;
    Legalize(tcx, t);
    t.ClearDelunayEdges();
    return ot;
}

Point& Sweep::NextFlipPoint(Point& ep, Point& eq, Triangle& ot, Point& op)
{
    Orientation o2d = Orient2d(eq, op, ep);
    if (o2d == CW) {
        // Right
        return *ot.PointCCW(op);
    }
    else if (o2d == CCW) {
        // Left
        return *ot.PointCW(op);
    }
    throw std::runtime_error("[Unsupported] Opposing point on constrained edge");
}

void Sweep::FlipScanEdgeEvent(SweepContext& tcx, Point& ep, Point& eq, Triangle& flip_triangle,
    Triangle& t, Point& p)
{
    Triangle& ot = t.NeighborAcross(p);
    Point& op = *ot.OppositePoint(t, p);

    if (InScanArea(eq, *flip_triangle.PointCCW(eq), *flip_triangle.PointCW(eq), op)) {
        // flip with new edge op->eq
        FlipEdgeEvent(tcx, eq, op, &ot, op);
        // TODO: Actually I just figured out that it should be possible to
        //       improve this by getting the next ot and op before the the above
        //       flip and continue the flipScanEdgeEvent here
        // set new ot and op here and loop back to inScanArea test
        // also need to set a new flip_triangle first
        // Turns out at first glance that this is somewhat complicated
        // so it will have to wait.
    }
    else{
        Point& newP = NextFlipPoint(ep, eq, ot, op);
        FlipScanEdgeEvent(tcx, ep, eq, flip_triangle, ot, newP);
    }
}

Sweep::~Sweep() {

    // Clean up memory
    for (size_t i = 0; i < nodes_.size(); i++) {
        delete nodes_[i];
    }

}

class CDT
{
public:

    /**
    * Constructor - add polyline with non repeating points
    *
    * @param polyline
    */
    CDT(const std::vector<Point*>& polyline);

    /**
    * Destructor - clean up memory
    */
    ~CDT();

    /**
    * Add a hole
    *
    * @param polyline
    */
    void AddHole(const std::vector<Point*>& polyline);

    /**
    * Add a steiner point
    *
    * @param point
    */
    void AddPoint(Point* point);

    /**
    * Triangulate - do this AFTER you've added the polyline, holes, and Steiner points
    */
    void Triangulate();

    /**
    * Get CDT triangles
    */
    std::vector<Triangle*> GetTriangles();

    /**
    * Get triangle map
    */
    std::list<Triangle*> GetMap();

private:

    /**
    * Internals
    */

    SweepContext* sweep_context_;
    Sweep* sweep_;

};

CDT::CDT(const std::vector<Point*>& polyline)
{
    sweep_context_ = new SweepContext(polyline);
    sweep_ = new Sweep;
}

void CDT::AddHole(const std::vector<Point*>& polyline)
{
    sweep_context_->AddHole(polyline);
}

void CDT::AddPoint(Point* point) {
    sweep_context_->AddPoint(point);
}

void CDT::Triangulate()
{
    sweep_->Triangulate(*sweep_context_);
}

std::vector<Triangle*> CDT::GetTriangles()
{
    return sweep_context_->GetTriangles();
}

std::list<Triangle*> CDT::GetMap()
{
    return sweep_context_->GetMap();
}

CDT::~CDT()
{
    delete sweep_context_;
    delete sweep_;
}


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

    vector<Integers> graph;
    graph.resize(m, vector<int>());
    vector<bool> used;
    used.resize(m);
    queue<int> qu;

    vector<Point> points(m);
    vector<Point*> ppoints(m);
    for (int i = 0; i < m; ++i) {
        points[i] = Point(dpoints[i]._x, dpoints[i]._y, i);
        ppoints[i] = &points[i];
    }
    CDT cdt(ppoints);

    std::vector<Triangle*> triangles;
    if (m > 2) {
        cdt.Triangulate();
        triangles = cdt.GetTriangles();
    }

    graph.resize(m);
    for (size_t i = 0; i < triangles.size(); ++i) {
        Triangle& t = *triangles[i];
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                if (j != k) {
                    if (t.GetPoint(j)->index >= 0 && t.GetPoint(k)->index >= 0) {
                        graph[t.GetPoint(j)->index].push_back(t.GetPoint(k)->index);
                    }
                }
            }
        }
    }

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
        TBase mind = INF;
        if (m > 10) {
            temp.clear();
            const Triangle& t = *(triangles[0]);
            for (int i = 0; i < 3; ++i) {
                if (t.GetPoint(i)->index >= 0) {
                    temp.push_back(t.GetPoint(0)->index);
                }
            }
            for (int k = 0; k < temp.size(); ++k) {
                qu.push(temp[k]);
                used[temp[k]] = true;
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
                    for (size_t k = 0; k < graph[index].size(); ++k) {
                        int v = graph[index][k];
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
        }
        else {
            dindices.resize(m);
            for (int j = 0; j < m; ++j) {
                dindices[j] = j;
            }
        }

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