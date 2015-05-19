#define _CRT_SECURE_NO_WARNINGS
#define M_PI 3.14159265358979323846

#include <cstdio>
#include <memory.h>

#include <vector>
#include <algorithm>

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

struct TPoint
{
	double x;
	double y;

	TPoint()
    {
    }

	TPoint(const TPoint& p)
        : x(p.x)
        , y(p.y) 
    {
    }

	TPoint(double _x, double _y) 
        : x(_x)
        , y(_y) 
    {
    
    }
};

struct TSite
{
	TPoint p;

	struct THalfEdge* edges;
};

struct TEdge
{
	TSite* 		sites[2];
	TPoint		pos[2];
	double		a;
	double		b;
	double		c;
	TEdge*		next;

	void create(struct TSite* s1, struct TSite* s2);
	void clipline(double width, double height);
};

class TVoronoi
{
public:
	TVoronoi();
	~TVoronoi();

	size_t get_required_mem();

	void generate(size_t num_sites, const TPoint* sites, int _width, int _height);

	const TSite* get_cells() const;
	const TEdge* get_edges() const;

private:
	void cleanup();

	TSite*	nextsite();
	void 	site_event(TSite* site);
	void 	circle_event();
	bool 	check_circle_event(THalfEdge* he1, THalfEdge* he2, TPoint& vertex);

	THalfEdge* get_edge_above_x(const TPoint& p);
	bool edge_intersect(const THalfEdge& e1, const THalfEdge& e2, TPoint& out) const;
	void endpos(struct TEdge* e, const TPoint& p, int direction);

	struct TEdge*		edges;
	struct THalfEdge*	beachline_start;
	struct THalfEdge*	beachline_end;
	struct THalfEdge*	last_inserted;
	struct TPriorityQueue* eventqueue;

	struct TSite*		sites;
	struct TSite*		bottomsite;
	int					numsites;
	int					currentsite;

	double 				width;
	double 				height;

	struct TEdge*		edgemem;
	struct TEdge*		edgepool;

	struct THalfEdge*	halfedgemem;
	struct THalfEdge*	halfedgepool;

	void**				eventmem;

	struct TEdge*		new_edge(TSite* s1, TSite* s2);
	void				delete_edge(TEdge*);

	struct THalfEdge*	new_halfedge(TEdge* e, int direction);
	void				delete_halfedge(THalfEdge*);
};

typedef bool (*FPriorityQueueCompare)(const void* node1, const void* node2);
typedef void (*FPriorityQueueSetpos)(const void* node, int pos);
typedef int  (*FPriorityQueueGetpos)(const void* node);
typedef int  (*FPriorityQueuePrint)(const void* node, int pos);

struct TPriorityQueue
{
	FPriorityQueueCompare	compare;
	FPriorityQueueSetpos	setpos;
	FPriorityQueueGetpos	getpos;
	int						maxnumitems;
	int						numitems;
	const void**			items;
};

void 		pq_create(TPriorityQueue* pq, int capacity, const void** buffer,
						FPriorityQueueCompare cmp,
						FPriorityQueueSetpos setpos,
						FPriorityQueueGetpos getpos);
bool 		pq_empty(TPriorityQueue* pq);
int 		pq_push(TPriorityQueue* pq, const void* node);
const void* pq_pop(TPriorityQueue* pq);
const void* pq_top(TPriorityQueue* pq);
void 		pq_remove(TPriorityQueue* pq, const void* node);
void		pq_print(TPriorityQueue* pq, FPriorityQueuePrint printfn);

static const int DIRECTION_LEFT  = 0;
static const int DIRECTION_RIGHT = 1;
static const double INVALID_VALUE = double(-1);

struct THalfEdge
{
	//Site*		site;
	TEdge*		edge;
	THalfEdge*	left;
	THalfEdge*	right;
	TPoint		vertex;
	double		y;
	int 		direction; // 0=left, 1=right
	int			pqpos;

	void			create(TEdge* edge, int direction);
	TSite*	        leftsite() const;
	TSite*	        rightsite() const;
	bool		 	rightof(const TPoint& p) const;

	void link(THalfEdge* newedge)
	{
		newedge->left = this;
		newedge->right = this->right;
		this->right->left = newedge;
		this->right = newedge;
	}

	void unlink()
	{
		left->right = right;
		right->left = left;
		left  = 0;
		right = 0;
	}
};

void THalfEdge::create(TEdge* e, int dir)
{
	edge 		= e;
	left 		= 0;
	right		= 0;
	direction 	= dir;
	pqpos		= 0;
	y			= 0;
}

TSite* THalfEdge::leftsite() const
{
	return edge->sites[direction];
}

TSite* THalfEdge::rightsite() const
{
	return edge ? edge->sites[1-direction] : 0;
}

bool THalfEdge::rightof(const TPoint& p) const
{
	const THalfEdge* he 		= this;
	const TEdge*		e		= he->edge;
	const TSite*		topsite = e->sites[1];

	bool right_of_site = p.x > topsite->p.x;
	if (right_of_site && he->direction == DIRECTION_LEFT)
		return true;
	if (!right_of_site && he->direction == DIRECTION_RIGHT)
		return false;

	double dxp, dyp, dxs, t1, t2, t3, yl;

	int above;
	if (e->a == double(1))
	{
		dyp = p.y - topsite->p.y;
		dxp = p.x - topsite->p.x;
		int fast = 0;
		if( (!right_of_site & (e->b < double(0))) | (right_of_site & (e->b >= double(0))) )
		{
			above = dyp >= e->b * dxp;
			fast = above;
		}
		else
		{
			above = p.x + p.y * e->b > e->c;
			if (e->b < double(0))
				above = !above;
			if (!above)
				fast = 1;
		}
		if (!fast)
		{
			dxs = topsite->p.x - e->sites[0]->p.x;
			above = e->b * (dxp * dxp - dyp * dyp)
					< dxs * dyp * (double(1) + double(2) * dxp / dxs + e->b * e->b);
			if (e->b < double(0))
				above = !above;
		}
	}
	else // e->b == 1
	{
		yl = e->c - e->a * p.x;
		t1 = p.y - yl;
		t2 = p.x - topsite->p.x;
		t3 = yl - topsite->p.y;
		above = t1 * t1 > t2 * t2 + t3 * t3;
	};
	return (he->direction == DIRECTION_LEFT ? above : !above);
}

void TEdge::create(TSite* s1, TSite* s2)
{
	TEdge* e = this;
	e->next = 0;
	sites[0] = s1;
	sites[1] = s2;
	pos[0].x = INVALID_VALUE;
	pos[1].x = INVALID_VALUE;

	// Create line equation between S1 and S2:
	// real_t a = -1 * (s2->p.y - s1->p.y);
	// real_t b = s2->p.x - s1->p.x;
	// //real_t c = -1 * (s2->p.x - s1->p.x) * s1->p.y + (s2->p.y - s1->p.y) * s1->p.x;
	//
	// // create perpendicular line
	// real_t pa = b;
	// real_t pb = -a;
	// //real_t pc = pa * s1->p.x + pb * s1->p.y;
	//
	// // Move to the mid point
	// real_t mx = s1->p.x + dx * real_t(0.5);
	// real_t my = s1->p.y + dy * real_t(0.5);
	// real_t pc = ( pa * mx + pb * my );

	double dx = s2->p.x - s1->p.x;
	double dy = s2->p.y - s1->p.y;

	// Simplify it, using dx and dy
	e->c = dx * (s1->p.x + dx * double(0.5)) + dy * (s1->p.y + dy * double(0.5));

	if (abs(dx) > abs(dy))
	{
		e->a = double(1);
		e->b = dy / dx;
		e->c /= dx;
	}
	else
	{
		e->a = dx / dy;
		e->b = double(1);
		e->c /= dy;
	}
}

void TEdge::clipline(double width, double height)
{
	TEdge* e = this;

	double pxmin = 0;
	double pxmax = width;
	double pymin = 0;
	double pymax = height;

	double x1, y1, x2, y2;
	TPoint* s1;
	TPoint* s2;
	if (e->a == double(1) && e->b >= double(0))
	{
		s1 = e->pos[1].x != INVALID_VALUE ? &e->pos[1] : 0;
		s2 = e->pos[0].x != INVALID_VALUE ? &e->pos[0] : 0;
	}
	else
	{
		s1 = e->pos[0].x != INVALID_VALUE ? &e->pos[0] : 0;
		s2 = e->pos[1].x != INVALID_VALUE ? &e->pos[1] : 0;
	};

	if (e->a == double(1))
	{
		y1 = pymin;
		if (s1 != 0 && s1->y > pymin)
		{
			y1 = s1->y;
		}
		if( y1 > pymax )
		{
			y1 = pymax;
		}
		x1 = e->c - e->b * y1;
		y2 = pymax;
		if (s2 != 0 && s2->y < pymax)
			y2 = s2->y;

		if( y2 < pymin )
		{
			y2 = pymin;
		}
		x2 = (e->c) - (e->b) * y2;
		if( ((x1 > pxmax) & (x2 > pxmax)) | ((x1 < pxmin) & (x2 < pxmin)) )
		{
			return;
		}
		if (x1 > pxmax)
		{
			x1 = pxmax;
			y1 = (e->c - x1) / e->b;
		}
		else if (x1 < pxmin)
		{
			x1 = pxmin;
			y1 = (e->c - x1) / e->b;
		}
		if (x2 > pxmax)
		{
			x2 = pxmax;
			y2 = (e->c - x2) / e->b;
		}
		else if (x2 < pxmin)
		{
			x2 = pxmin;
			y2 = (e->c - x2) / e->b;
		}
	}
	else
	{
		x1 = pxmin;
		if( s1 != 0 && s1->x > pxmin )
			x1 = s1->x;
		if( x1 > pxmax )
		{
			x1 = pxmax;
		}
		y1 = e->c - e->a * x1;
		x2 = pxmax;
		if( s2 != 0 && s2->x < pxmax )
			x2 = s2->x;
		if( x2 < pxmin )
		{
			x2 = pxmin;
		}
		y2 = e->c - e->a * x2;
		if( ((y1 > pymax) & (y2 > pymax)) | ((y1 < pymin) & (y2 < pymin)) )
		{
			return;
		}
		if( y1 > pymax )
		{
			y1 = pymax;
			x1 = (e->c - y1) / e->a;
		}
		else if( y1 < pymin )
		{
			y1 = pymin;
			x1 = (e->c - y1) / e->a;
		}
		if( y2 > pymax )
		{
			y2 = pymax;
			x2 = (e->c - y2) / e->a;
		}
		else if( y2 < pymin )
		{
			y2 = pymin;
			x2 = (e->c - y2) / e->a;
		}
	};

	pos[0].x = x1;
	pos[0].y = y1;
	pos[1].x = x2;
	pos[1].y = y2;
}

static int pq_moveup(TPriorityQueue* pq, int pos)
{
	const void* node = pq->items[pos];

	for( int parent = (pos >> 1);
		 pos > 1 && pq->compare(pq->items[parent], node);
		 pos = parent, parent = parent >> 1)
	{
		pq->items[pos] = pq->items[parent];
		pq->setpos( (void*)pq->items[pos], pos );
	}

	pq->items[pos] = node;
	pq->setpos( (void*)pq->items[pos], pos );
	return pos;
}

static int pq_maxchild(TPriorityQueue* pq, int pos)
{
	int child = pos << 1;
	if( child >= pq->numitems )
		return 0;
	if( (child + 1) < pq->numitems && pq->compare(pq->items[child], pq->items[child+1]) )
		return child+1;
	return child;
}

static int pq_movedown(TPriorityQueue* pq, int pos)
{
	const void* node = pq->items[pos];

	int child;
	while( (child = pq_maxchild(pq, pos)) &&
			pq->compare( node, pq->items[child] ) )
	{
		pq->items[pos] = pq->items[child];
		pq->setpos( (void*)pq->items[pos], pos );
		pos = child;
	}

	pq->items[pos] = node;
	pq->setpos( (void*)pq->items[pos], pos );
	return pos;
}

void pq_create(TPriorityQueue* pq, int capacity, const void** buffer,
				FPriorityQueueCompare cmp,
				FPriorityQueueSetpos setpos,
				FPriorityQueueGetpos getpos)
{
	pq->compare 	= cmp;
	pq->setpos		= setpos;
	pq->getpos		= getpos;
	pq->maxnumitems = capacity;
	pq->numitems	= 1;
	pq->items 		= buffer;
}

bool pq_empty(TPriorityQueue* pq)
{
	return pq->numitems == 1;
}

int pq_push(TPriorityQueue* pq, const void* node)
{
	int n = pq->numitems++;
	pq->items[n] = node;
	return pq_moveup(pq, n);
}

const void* pq_pop(TPriorityQueue* pq)
{
	if( pq->numitems == 1 )
		return 0;

	const void* node = pq->items[1];
	pq->items[1] = pq->items[--pq->numitems];
	pq_movedown(pq, 1);
	return node;
}

const void* pq_top(TPriorityQueue* pq)
{
	if( pq->numitems == 1 )
		return 0;
	return pq->items[1];
}

void pq_remove(TPriorityQueue* pq, const void* node)
{
	if( pq->numitems == 1 )
		return;
	int pos = pq->getpos(node);
	if( pos == 0 )
		return;

	pos = 1;
	for( ; pos < pq->numitems; ++pos )
	{
		if( pq->items[pos] == node )
			break;
	}

	pq->items[pos] = pq->items[--pq->numitems];
	if( pq->compare( node, pq->items[pos] ) )
		pq_moveup( pq, pos );
	else
		pq_movedown( pq, pos );
	pq->setpos( (void*)node, 0 );
}

void pq_print(TPriorityQueue* pq, FPriorityQueuePrint printfn)
{
	printf("\tPQ\n");
	for( int i = 1; i < pq->numitems; ++i )
	{
		printfn(pq->items[i], i);
	}
	printf("\t-----\n");
}

static void he_print(const THalfEdge* he, int pos)
{
	printf("\t%g, %g  y: %g\n", he->vertex.x, he->vertex.y, he->y);
}

static inline bool he_compare(const THalfEdge* he1, const THalfEdge* he2)
{
	return (he1->y > he2->y) || ((he1->y == he2->y) && (he1->vertex.x > he2->vertex.x));
}

static inline void he_setpos(THalfEdge* he, int pos)
{
	he->pqpos = pos;
}

static inline int he_getpos(const THalfEdge* he)
{
	return he->pqpos;
}

static inline int point_cmp(const void *p1, const void *p2)
{
	const TPoint& s1 = *(TPoint*)p1;
	const TPoint& s2 = *(TPoint*)p2;
	if (s1.y != s2.y)
    {
        if (s1.y > s2.y)
            return 1;
        else
            return -1;
    }
    if (s1.x > s2.x)
        return 1;
    else if (s1.x < s2.x)
        return -1;
    return 0;
}

static inline bool pt_less(const TPoint& pt1, const TPoint& pt2)
{
	return (pt1.y < pt2.y) || ((pt1.y == pt2.y) && (pt1.x < pt2.x));
}

static inline double pt_dist(const TPoint& pt1, const TPoint& pt2)
{
	double dx = pt1.x - pt2.x;
	double dy = pt1.y - pt2.y;
	return (double) sqrt(dx*dx + dy*dy);
}

static void printbeach(const THalfEdge* he)
{
	printf("(%g, %g, y: %g), ", he->vertex.x, he->vertex.y, he->y );
	if( he->right )
		printbeach(he->right);
}

TVoronoi::TVoronoi()
{
	sites 			= 0;
	beachline_start = 0;
	beachline_end	= 0;
	edges			= 0;
	eventmem		= 0;
	eventqueue		= 0;
	edgemem			= 0;
	halfedgemem		= 0;
}


void TVoronoi::cleanup()
{
	if( sites )
		delete[] sites;
	if( eventmem )
		free(eventmem);
	if( edgemem )
		delete[] edgemem;
	if( halfedgemem )
		delete[] halfedgemem;
	if( eventqueue )
		free(eventqueue);
	sites 			= 0;
	beachline_start = 0;
	beachline_end	= 0;
	edges			= 0;
	eventqueue		= 0;
	eventmem		= 0;
}

TVoronoi::~TVoronoi()
{
	cleanup();
}

TSite* TVoronoi::nextsite()
{
	return (currentsite < numsites) ? &sites[currentsite++] : 0;
}

void TVoronoi::generate(size_t num_sites, const TPoint* _sites, int _width, int _height)
{
	cleanup();

	sites = new TSite[num_sites];

	size_t max_num_edges = num_sites * 3 - 6 + 1;

	edgemem		= new TEdge[max_num_edges];
	edgepool	= edgemem;
	for (int i = 0; i < max_num_edges-1; ++i)
	{
		edgemem[i].next = &edgemem[i+1];
	}
	edgemem[max_num_edges-1].next = 0;

	size_t max_num_arcs		 = 2 * num_sites - 1;
	size_t max_num_halfedges = max_num_arcs * 2 + 2; // Each arc gets 2 half edges, and we want 2 extra at the start/end

	halfedgemem 	= new THalfEdge[max_num_halfedges];
	halfedgepool 	= halfedgemem;

	memset(halfedgemem, 0, max_num_halfedges * sizeof(THalfEdge));

	for (int i = 0; i < max_num_halfedges-1; ++i)
	{
		halfedgemem[i].right = &halfedgemem[i+1];
	}
	halfedgemem[max_num_halfedges-1].right = 0;

	beachline_start = halfedgepool;
	halfedgepool = halfedgepool->right;
	beachline_end = halfedgepool;
	halfedgepool = halfedgepool->right;

	beachline_start->left 	= 0;
	beachline_start->right 	= beachline_end;
	beachline_end->left		= beachline_start;
	beachline_end->right	= 0;

	last_inserted = 0;

	int max_num_events = num_sites * 2 + 1;
	eventmem = (void**)malloc(max_num_events * sizeof(void*));

	eventqueue = (TPriorityQueue*)malloc(sizeof(TPriorityQueue));
	pq_create(eventqueue, max_num_events, (const void**)eventmem,
				(FPriorityQueueCompare)he_compare,
				(FPriorityQueueSetpos)he_setpos,
				(FPriorityQueueGetpos)he_getpos);


	for( size_t i = 0; i < num_sites; ++i )
	{
		sites[i].p 	= _sites[i];
		//sites[i].id = (int)i;
        sites[i].edges = 0;
	}

	// Remove duplicates, to avoid anomalies
	qsort(sites, num_sites, sizeof(TSite), point_cmp);

	unsigned int offset = 0;
	for (int is = 1; is < num_sites; is++)
	{
		if( sites[is].p.y == sites[is - 1].p.y && sites[is].p.x == sites[is - 1].p.x )
		{
			offset++;
			continue;
		}
		else if (offset > 0)
		{
			sites[is - offset] = sites[is];
		}
	}
	num_sites 	-= offset;
	numsites 	= num_sites;
	currentsite = 0;

	width = _width;
	height = _height;

	bottomsite = nextsite();

	TSite* site = nextsite();

	while (true)
	{
		TPoint lowest_pq_point;
		if( !pq_empty(eventqueue) )
		{
			THalfEdge* he = (THalfEdge*)pq_top(eventqueue);
			lowest_pq_point.x = he->vertex.x;
			lowest_pq_point.y = he->y;
		}

		if( site != 0 && (pq_empty(eventqueue) || pt_less(site->p, lowest_pq_point) ) )
		{
			site_event(site);
			site = nextsite();
		}
		else if( !pq_empty(eventqueue) )
		{
			circle_event();
		}
		else
		{
			break;
		}

		//pq_print(eventqueue, (FPriorityQueuePrint)he_print);
	}

	TEdge* edge = edges;
	while (edge)
	{
		edge->clipline(width, height);
		edge = edge->next;
	}
}

const TSite* TVoronoi::get_cells() const
{
	return sites;
}

const TEdge* TVoronoi::get_edges() const
{
	return edges;
}

void TVoronoi::site_event(TSite* site)
{
	//printf("%s   %g, %g\n", __FUNCTION__, site->p.x, site->p.y);

	THalfEdge* left 	 = get_edge_above_x(site->p);
	THalfEdge* right	 = left->right;
	TSite*	  bottom     = left->rightsite();
	if (!bottom)
		bottom = bottomsite;

	TEdge* edge = new_edge(bottom, site);
	edge->next = edges;
	edges = edge;

	THalfEdge* edge1 = new_halfedge(edge, DIRECTION_LEFT);
	THalfEdge* edge2 = new_halfedge(edge, DIRECTION_RIGHT);

	THalfEdge* leftleft = left->left;

	left->link(edge1);
	edge1->link(edge2);

	last_inserted = edge1;

	TPoint p;
	if (check_circle_event(left, edge1, p))
	{
		if( left->pqpos )
			pq_remove(eventqueue, left);
		left->vertex 	= p;
		left->y		 	= p.y + pt_dist(site->p, p);
		pq_push(eventqueue, left);
	}
	if (check_circle_event(edge2, right, p))
	{
		edge2->vertex	= p;
		edge2->y		= p.y + pt_dist(site->p, p);
		pq_push(eventqueue, edge2);
	}
}

void TVoronoi::circle_event()
{
	THalfEdge* left = (THalfEdge*)pq_pop(eventqueue);

	//printf("%s   %g, %g  Y: %g\n", __FUNCTION__, left->vertex.x, left->vertex.y, left->y);
	//pq_print(eventqueue, (FPriorityQueuePrint)he_print);

	THalfEdge* leftleft 	= left->left;
	THalfEdge* right		= left->right;
	THalfEdge* rightright   = right->right;
	TSite* bottom = left->leftsite();
	TSite* top 	 = right->rightsite();

	TPoint vertex = left->vertex;
	endpos(left->edge, vertex, left->direction);
	endpos(right->edge, vertex, right->direction);

	if( last_inserted == left )
		last_inserted = leftleft;
	else if( last_inserted == right )
		last_inserted = rightright;
	pq_remove(eventqueue, right);
	left->unlink();
	right->unlink();
    delete_halfedge(left);
    delete_halfedge(right);

	int direction = DIRECTION_LEFT;
	if (bottom->p.y > top->p.y)
	{
		TSite* temp = bottom;
		bottom = top;
		top = temp;
		direction = DIRECTION_RIGHT;
	}

	TEdge* edge = new_edge(bottom, top);
	edge->next = edges;
	edges = edge;

	THalfEdge* he = new_halfedge(edge, direction);
	leftleft->link(he);
	endpos(edge, vertex, DIRECTION_RIGHT - direction);

	TPoint p;
	if (check_circle_event(leftleft, he, p))
	{
		if( leftleft->pqpos )
			pq_remove(eventqueue, leftleft);
		leftleft->vertex 	= p;
		leftleft->y		 	= p.y + pt_dist(bottom->p, p);
		pq_push(eventqueue, leftleft);
	}
	if (check_circle_event(he, rightright, p))
	{
		he->vertex 		= p;
		he->y		 	= p.y + pt_dist(bottom->p, p);
		pq_push(eventqueue, he);
	}
}

void TVoronoi::endpos(TEdge* e, const TPoint& p, int direction)
{
	e->pos[direction] = p;

	if( e->pos[0].x != double(-1) && e->pos[1].x != double(-1) )
	{
		e->clipline(width, height);

		// TODO: Add to site

		//e->next = edges;
		//edges = e;
	}
}

bool TVoronoi::edge_intersect(const THalfEdge& he1, const THalfEdge& he2, TPoint& out) const
{
	const TEdge& e1 = *he1.edge;
	const TEdge& e2 = *he2.edge;

	double dx = e2.sites[1]->p.x - e1.sites[1]->p.x;
	double dy = e2.sites[1]->p.y - e1.sites[1]->p.y;

	if( dx == 0 && dy == 0 )
	{
		return false;
	}

	double d = e1.a * e2.b - e1.b * e2.a;
	if(abs(d) < double(0.00001))
	{
		return false;
	}
	out.x = (e1.c * e2.b - e1.b * e2.c) / d;
	out.y = (e1.a * e2.c - e1.c * e2.a) / d;

	const TEdge* e;
	const THalfEdge* he;
	if (pt_less( e1.sites[1]->p, e2.sites[1]->p))
	{
		he = &he1;
		e = &e1;
	}
	else
	{
		he = &he2;
		e = &e2;
	}

	int right_of_site = out.x >= e->sites[1]->p.x;
	if ((right_of_site && he->direction == DIRECTION_LEFT) || (!right_of_site && he->direction == DIRECTION_RIGHT))
	{
		return false;
	}

	return true;
}

bool TVoronoi::check_circle_event(THalfEdge* he1, THalfEdge* he2, TPoint& vertex)
{
	TEdge* e1 = he1->edge;
	TEdge* e2 = he2->edge;
	if( e1 == 0 || e2 == 0 || e1->sites[1] == e2->sites[1] )
	{
		return false;
	}

	return edge_intersect(*he1, *he2, vertex);
}

struct THalfEdge* TVoronoi::get_edge_above_x(const TPoint& p)
{
	// Gets the arc on the beach line at the x coordinate (i.e. right above the new site event)

	// A good guess it's close by (Can be optimized)
	THalfEdge* he = last_inserted;
    if (!he)
    {
        if( p.x < width / 2 )
            he = beachline_start;
        else
            he = beachline_end;
    }

	//
	if( he == beachline_start || (he != beachline_end && he->rightof(p)) )
	{
		do {
			
            he = he->right;
		}
		while( he != beachline_end && he->rightof(p) );

		he = he->left;
	}
	else
	{
		do
        {
			he = he->left;
		}
		while( he != beachline_start && !he->rightof(p) );
	}

	return he;
}

struct TEdge* TVoronoi::new_edge(TSite* s1, TSite* s2)
{
	TEdge* p = edgepool;
	edgepool = edgepool->next;
	p->create(s1, s2);
	return p;
}

void TVoronoi::delete_edge(TEdge* e)
{
	e->next = edgepool;
	edgepool = e;
}

THalfEdge* TVoronoi::new_halfedge(TEdge* e, int direction)
{
	THalfEdge* he = halfedgepool;
	halfedgepool = halfedgepool->right;
	he->create(e, direction);
	return he;
}

void TVoronoi::delete_halfedge(THalfEdge* he)
{
	he->right = halfedgepool;
    halfedgepool = he;
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

void GenRandom()
{
    FILE* fOut = fopen("random.txt", "w");
    static const int N = 100000;
    fprintf(fOut, "%d\n", N);
    for (int i = 0; i < N; ++i)
    {
        fprintf(fOut, "%lf %lf\n", (rand() % 1000000)/1000000, (rand() % 1000000)/1000000);
    }
    static const int M = 10000;
    fprintf(fOut, "%d\n", M);
    for (int i = 0; i < M; ++i)
    {
        fprintf(fOut, "%lf %lf\n", (rand() % 1000000)/1000000, (rand() % 1000000)/1000000);
    }
    fclose(fOut);
}

int main()
{
#ifndef ONLINE_JUDGE
    // GenBig();
    // freopen("big.txt", "r", stdin);
    GenRandom();
    freopen("random.txt", "r", stdin);
    // freopen("input.txt", "r", stdin);
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

    typedef vector<TPoint> TPoints;
    TPoints points(m);
    for (int i = 0; i < m; ++i)
    {
        points[i] = TPoint(xc[i], yc[i]);
    }

    TVoronoi v;
    v.generate(points.size(), &points[0], 100, 100);
   
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

        result.clear();
        for (int i = 0; i < m; ++i)
        {
            result.push_back(i);
        }

        result2.clear();
        long double mind = 1e15;
        long double mindMin = mind;
        long double mindMax = mind;
        for (auto realindex : result)
        {
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

        Output(result2);
    }

    return 0;
}