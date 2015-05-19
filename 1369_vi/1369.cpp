#define _CRT_SECURE_NO_WARNINGS

#include <immintrin.h>
#include <intrin.h>

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

struct Triad
{
    int a, b, c;
    int ab, bc, ac;  // adjacent edges index to neighbouring triangle.
    TBase ro, R, C;
    
    Triad() {
    }
    
    Triad(int x, int y) 
        : a(x)
        , b(y)
        , c(0)
        , ab(-1)
        , bc(-1)
        , ac(-1)
        , ro(-1)
        , R(0)
        , C(0)
    {
    }

    Triad(int x, int y, int z) 
        : a(x)
        , b(y)
        , c(z)
        , ab(-1)
        , bc(-1)
        , ac(-1)
        , ro(-1)
        , R(0)
        , C(0)
    {
    }

    Triad(const Triad& p) 
        : a(p.a)
        , b(p.b)
        , c(p.c)
        , ab(p.ab)
        , bc(p.bc)
        , ac(p.ac)
        , ro(p.ro)
        , R(p.R)
        , C(p.C)
    {
    }

    Triad& operator=(const Triad& p) {
        a = p.a;
        b = p.b;
        c = p.c;

        ab = p.ab;
        bc = p.bc;
        ac = p.ac;

        ro = p.ro;
        R = p.R;
        C = p.C;

        return *this;
    }
};


struct Shx
{
    int id, trid;
    TBase r, c, tr, tc;
    TBase ro;
    
    Shx() {
    }

    Shx(TBase a, TBase b)
        : r(a)
        , c(b)
        , ro(0.0)
        , tr(0.0)
        , tc(0.0)
        , id(-1)
    {
    }

    Shx(TBase a, TBase b, TBase x)
        : r(a)
        , c(b)
        , ro(x)
        , id(-1)
        , tr(0)
        , tc(0)
    {
    }

    Shx(const Shx& p) 
        : id(p.id)
        , trid(p.trid)
        , r(p.r)
        , c(p.c)
        , tr(p.tr)
        , tc(p.tc)
        , ro(p.ro)
    {
    }

    Shx& operator=(const Shx& p) {
        id = p.id;
        trid = p.trid;
        r = p.r;
        c = p.c;
        tr = p.tr;
        tc = p.tc;
        ro = p.ro;
        return *this;
    }
};

inline bool operator<(const Shx& a, const Shx& b) {
    if (a.ro == b.ro) {
        if (a.r == b.r) {
            return a.c < b.c;
        }
        return a.r < b.r;
    }
    return a.ro <  b.ro;
}

struct Dupex
{
    int id;
    TBase r, c;

    Dupex() {
    }

    Dupex(TBase a, TBase b)
        : r(a)
        , c(b)
        , id(-1)
    {
    }

    Dupex(TBase a, TBase b, int x)
        : r(a)
        , c(b)
        , id(x)
    {
    }

    Dupex(const Dupex& p)
        : id(p.id)
        , r(p.r)
        , c(p.c)
    {
    }

    Dupex& operator=(const Dupex& p) {
        id = p.id;
        r = p.r;
        c = p.c;
        return *this;
    }
};

// sort into descending order (for use in corner responce ranking).
inline bool operator<(const Dupex& a, const Dupex& b) {
    if (a.r == b.r) {
        return a.c < b.c;
    }
    return a.r < b.r;
}


int s_hull_pro(std::vector<Shx>& pts, std::vector<Triad>& triads);
void circle_cent2(TBase r1, TBase c1, TBase r2, TBase c2, TBase r3, TBase c3, TBase& r, TBase& c, TBase& ro2);
void circle_cent4(TBase r1, TBase c1, TBase r2, TBase c2, TBase r3, TBase c3, TBase& r, TBase& c, TBase& ro2);
int Cline_Renka_test(TBase& Ax, TBase& Ay, TBase& Bx, TBase& By, TBase& Cx, TBase& Cy, TBase& Dx, TBase& Dy);
int T_flip_pro(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, int numt, int start, std::set<int>& ids);
int T_flip_pro_idx(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, std::set<int>& ids);
int T_flip_edge(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, int numt, int start, std::set<int>& ids);
int  test_center(Shx& pt0, Shx& pt1, Shx& pt2);

void circle_cent2(TBase r1, TBase c1, TBase r2, TBase c2, TBase r3, TBase c3, TBase& r, TBase& c, TBase& ro2) {
    /*
    *  function to return the center of a circle and its radius
    * degenerate case should never be passed to this routine!!!!!!!!!!!!!
    * but will return r0 = -1 if it is.
    */

    TBase a1 = (r1 + r2) / 2.0;
    TBase a2 = (c1 + c2) / 2.0;
    TBase b1 = (r3 + r2) / 2.0;
    TBase b2 = (c3 + c2) / 2.0;

    TBase e2 = r1 - r2;
    TBase e1 = -c1 + c2;

    TBase q2 = r3 - r2;
    TBase q1 = -c3 + c2;

    r = 0; 
    c = 0; 
    ro2 = -1;
    if (e1*-q2 + e2*q1 == 0) 
        return;

    TBase beta = (-e2*(b1 - a1) + e1*(b2 - a2)) / (e2*q1 - e1*q2);

    r = b1 + q1*beta;
    c = b2 + q2*beta;

    ro2 = (r1 - r)*(r1 - r) + (c1 - c)*(c1 - c);    
}

/*  version in which the ids of the triangles associated with the sides of the hull are tracked.


*/

int s_hull_pro(std::vector<Shx>& pts, std::vector<Triad>& triads) {
    int nump = (int)pts.size();

    if (nump < 3) {
        throw Exception("less than 3 points, aborting ");
    }

    TBase r = pts[0].r;
    TBase c = pts[0].c;
    for (int k = 0; k < nump; ++k) {
        TBase dr = pts[k].r - r;
        TBase dc = pts[k].c - c;

        pts[k].ro = dr*dr + dc*dc;
    }

    sort(pts.begin(), pts.end());

    TBase r1 = pts[0].r;
    TBase c1 = pts[0].c;

    TBase r2 = pts[1].r;
    TBase c2 = pts[1].c;
    int mid = -1;
    TBase romin2 = 9.0e20, ro2, R, C;

    Shx pt0 = pts[0];
    Shx pt1 = pts[1];
    int k = 2;
    while (k < nump) {
        circle_cent2(r1, c1, r2, c2, pts[k].r, pts[k].c, r, c, ro2);
        if ((ro2 < romin2) && ro2 > 0 && (test_center(pt0, pt1, pts[k]) > 0)) {
            mid = k;
            romin2 = ro2;
            R = r;
            C = c;
        }
        else if (romin2 * 4 < pts[k].ro) {
            k = nump;
        }

        ++k;
    }

    if (mid < 0) {
        throw Exception("linear structure, aborting");
    }

    Shx pt2 = pts[mid];

    int ptest = test_center(pt0, pt1, pt2);
    if (ptest < 0) {
        throw Exception("warning: obtuce seed triangle selected");
    }

    pts.erase(pts.begin() + mid);  // necessary for round off reasons:((((((
    pts.erase(pts.begin());
    pts.erase(pts.begin());

    for (int k = 0; k < nump - 3; ++k) {
        TBase dr = pts[k].r - R;
        TBase dc = pts[k].c - C;

        pts[k].ro = dr*dr + dc*dc;
    }

    sort(pts.begin(), pts.end());

    pts.insert(pts.begin(), pt2);
    pts.insert(pts.begin(), pt1);
    pts.insert(pts.begin(), pt0);

    std::vector<int> slump;
    slump.resize(nump);

    for (int k = 0; k < nump; ++k) {
        if (pts[k].id < nump) {
            slump[pts[k].id] = k;
        } else {
            int mx = pts[k].id + 1;
            while ((int)slump.size() <= mx) {
                slump.push_back(0);
            }
            slump[pts[k].id] = k;
        }
    }

    std::vector<Shx> hull;

    r = (pts[0].r + pts[1].r + pts[2].r) / (float) 3.0;
    c = (pts[0].c + pts[1].c + pts[2].c) / (float) 3.0;

    TBase dr0 = pts[0].r - r;
    TBase dc0 = pts[0].c - c;
    TBase tr01 = pts[1].r - pts[0].r;
    TBase tc01 = pts[1].c - pts[0].c;

    TBase df = -tr01* dc0 + tc01*dr0;
    if (df < 0) {   // [ 0 1 2 ]
        pt0.tr = pt1.r - pt0.r;
        pt0.tc = pt1.c - pt0.c;
        pt0.trid = 0;
        hull.push_back(pt0);

        pt1.tr = pt2.r - pt1.r;
        pt1.tc = pt2.c - pt1.c;
        pt1.trid = 0;
        hull.push_back(pt1);

        pt2.tr = pt0.r - pt2.r;
        pt2.tc = pt0.c - pt2.c;
        pt2.trid = 0;
        hull.push_back(pt2);

        Triad tri(pt0.id, pt1.id, pt2.id);
        tri.ro = romin2;
        tri.R = R;
        tri.C = C;

        triads.push_back(tri);
    } else {          // [ 0 2 1 ] as anti-clockwise turning is the work of the devil....
        pt0.tr = pt2.r - pt0.r;
        pt0.tc = pt2.c - pt0.c;
        pt0.trid = 0;
        hull.push_back(pt0);

        pt2.tr = pt1.r - pt2.r;
        pt2.tc = pt1.c - pt2.c;
        pt2.trid = 0;
        hull.push_back(pt2);

        pt1.tr = pt0.r - pt1.r;
        pt1.tc = pt0.c - pt1.c;
        pt1.trid = 0;
        hull.push_back(pt1);

        Triad tri(pt0.id, pt2.id, pt1.id);
        tri.ro = romin2;
        tri.R = R;
        tri.C = C;
        
        triads.push_back(tri);
    }

    // add new points into hull (removing obscured ones from the chain)
    // and creating triangles....
    // that will need to be flipped.

    TBase dr, dc, rx, cx;
    Shx ptx;
    int numt;

    for (int k = 3; k < nump; ++k) {
        rx = pts[k].r;    
        cx = pts[k].c;
        ptx.r = rx;
        ptx.c = cx;
        ptx.id = pts[k].id;

        int numh = (int)hull.size();
        int numh_old = numh;
        dr = rx - hull[0].r;    
        dc = cx - hull[0].c;  // outwards pointing from hull[0] to pt.

        std::vector<int> pidx, tridx;
        int hidx;  // new hull point location within hull.....

        TBase df = -dc* hull[0].tr + dr*hull[0].tc;    // visibility test vector.
        if (df < 0) {  // starting with a visible hull facet !!!
            int e1 = 1;
            int e2 = numh;
            hidx = 0;

            // check to see if segment numh is also visible
            df = -dc* hull[numh - 1].tr + dr*hull[numh - 1].tc;
            //cerr << df << ' ' ;
            if (df < 0) {    // visible.
                pidx.push_back(hull[numh - 1].id);
                tridx.push_back(hull[numh - 1].trid);

                for (int h = 0; h < numh - 1; ++h) {
                    // if segment h is visible delete h
                    dr = rx - hull[h].r;    dc = cx - hull[h].c;
                    df = -dc* hull[h].tr + dr*hull[h].tc;
                    pidx.push_back(hull[h].id);
                    tridx.push_back(hull[h].trid);
                    if (df < 0) {
                        hull.erase(hull.begin() + h);
                        h--;
                        numh--;
                    } else {	  // quit on invisibility
                        ptx.tr = hull[h].r - ptx.r;
                        ptx.tc = hull[h].c - ptx.c;

                        hull.insert(hull.begin(), ptx);
                        numh++;
                        break;
                    }
                }
                // look backwards through the hull structure.

                for (int h = numh - 2; h>0; --h) {
                    // if segment h is visible delete h + 1
                    dr = rx - hull[h].r;    dc = cx - hull[h].c;
                    df = -dc* hull[h].tr + dr*hull[h].tc;

                    if (df < 0) {  // h is visible 
                        pidx.insert(pidx.begin(), hull[h].id);
                        tridx.insert(tridx.begin(), hull[h].trid);
                        hull.erase(hull.begin() + h + 1);  // erase end of chain
                    } else {
                        h = (int)hull.size() - 1;
                        hull[h].tr = -hull[h].r + ptx.r;   // points at start of chain.
                        hull[h].tc = -hull[h].c + ptx.c;
                        break;
                    }
                }

                df = 9;
            } else {
                //	cerr << df << ' ' << endl;
                hidx = 1;  // keep pt hull[0]
                tridx.push_back(hull[0].trid);
                pidx.push_back(hull[0].id);

                for (int h = 1; h<numh; h++) {
                    // if segment h is visible delete h  
                    dr = rx - hull[h].r;    dc = cx - hull[h].c;
                    df = -dc* hull[h].tr + dr*hull[h].tc;
                    pidx.push_back(hull[h].id);
                    tridx.push_back(hull[h].trid);
                    if (df < 0) {                     // visible
                        hull.erase(hull.begin() + h);
                        h--;
                        numh--;
                    } else {	  // quit on invisibility
                        ptx.tr = hull[h].r - ptx.r;
                        ptx.tc = hull[h].c - ptx.c;

                        hull[h - 1].tr = ptx.r - hull[h - 1].r;
                        hull[h - 1].tc = ptx.c - hull[h - 1].c;

                        hull.insert(hull.begin() + h, ptx);
                        break;
                    }
                }
            }

            df = 8;
        } else {
            int e1 = -1, e2 = numh;
            for (int h = 1; h<numh; h++) {
                dr = rx - hull[h].r;    
                dc = cx - hull[h].c;
                df = -dc* hull[h].tr + dr*hull[h].tc;
                if (df < 0) {
                    if (e1 < 0) 
                        e1 = h;  // fist visible
                } else {
                    if (e1 > 0) { // first invisible segment.
                        e2 = h;
                        break;
                    }
                }
            }


            // triangle pidx starts at e1 and ends at e2 (inclusive).	
            if (e2 < numh) {
                for (int e = e1; e <= e2; e++) {
                    pidx.push_back(hull[e].id);
                    tridx.push_back(hull[e].trid);
                }
            } else {
                for (int e = e1; e<e2; e++) {
                    pidx.push_back(hull[e].id);
                    tridx.push_back(hull[e].trid);   // there are only n-1 triangles from n hull pts.
                }
                pidx.push_back(hull[0].id);
            }


            // erase elements e1+1 : e2-1 inclusive.

            if (e1 < e2 - 1) {
                hull.erase(hull.begin() + e1 + 1, hull.begin() + e2);
            }

            // insert ptx at location e1+1.
            if (e2 == numh) {
                ptx.tr = hull[0].r - ptx.r;
                ptx.tc = hull[0].c - ptx.c;
            } else {
                ptx.tr = hull[e1 + 1].r - ptx.r;
                ptx.tc = hull[e1 + 1].c - ptx.c;
            }

            hull[e1].tr = ptx.r - hull[e1].r;
            hull[e1].tc = ptx.c - hull[e1].c;

            hull.insert(hull.begin() + e1 + 1, ptx);
            hidx = e1 + 1;
        }

        int a = ptx.id;
        int T0;
        Triad trx(a, 0, 0);
        r1 = pts[slump[a]].r;
        c1 = pts[slump[a]].c;

        int npx = (int)pidx.size() - 1;
        numt = (int)triads.size();
        T0 = numt;

        if (npx == 1) {
            trx.b = pidx[0];
            trx.c = pidx[1];

            trx.bc = tridx[0];
            trx.ab = -1;
            trx.ac = -1;

            // index back into the triads.
            Triad &txx = triads[tridx[0]];
            if ((trx.b == txx.a && trx.c == txx.b) || (trx.b == txx.b && trx.c == txx.a)) {
                txx.ab = numt;
            }
            else if ((trx.b == txx.a && trx.c == txx.c) || (trx.b == txx.c && trx.c == txx.a)) {
                txx.ac = numt;
            }
            else if ((trx.b == txx.b && trx.c == txx.c) || (trx.b == txx.c && trx.c == txx.b)) {
                txx.bc = numt;
            }


            hull[hidx].trid = numt;
            if (hidx > 0) {
                hull[hidx - 1].trid = numt;
            } else {
                numh = (int)hull.size();
                hull[numh - 1].trid = numt;
            }
            triads.push_back(trx);
            numt++;
        } else {
            trx.ab = -1;
            for (int p = 0; p < npx; ++p) {
                trx.b = pidx[p];
                trx.c = pidx[p + 1];

                trx.bc = tridx[p];
                if (p > 0)
                    trx.ab = numt - 1;
                trx.ac = numt + 1;

                // index back into the triads.
                Triad& txx = triads[tridx[p]];
                if ((trx.b == txx.a && trx.c == txx.b) || (trx.b == txx.b && trx.c == txx.a)) {
                    txx.ab = numt;
                }
                else if ((trx.b == txx.a && trx.c == txx.c) || (trx.b == txx.c && trx.c == txx.a)) {
                    txx.ac = numt;
                }
                else if ((trx.b == txx.b && trx.c == txx.c) || (trx.b == txx.c && trx.c == txx.b)) {
                    txx.bc = numt;
                }

                triads.push_back(trx);
                numt++;
            }
            triads[numt - 1].ac = -1;

            hull[hidx].trid = numt - 1;
            if (hidx > 0) {
                hull[hidx - 1].trid = T0;
            } else {
                numh = (int)hull.size();
                hull[numh - 1].trid = T0;
            }
        }
    }

    std::set<int> ids;
    std::set<int> ids2;

    int tf = T_flip_pro(pts, triads, slump, numt, 0, ids);
    if (tf < 0) {
        throw Exception("cannot triangualte this set ");

        return -3;
    }

    int nits = (int)ids.size();
    int nit = 1;
    while (nits > 0 && nit < 50) {
        tf = T_flip_pro_idx(pts, triads, slump, ids);
        nits = (int)ids.size();
        nit++;
        if (tf < 0) {
            throw Exception("cannot triangualte this set ");

            return -4;
        }
    }

    ids.clear();
    nits = T_flip_edge(pts, triads, slump, numt, 0, ids);
    nit = 0;

    while (nits > 0 && nit < 100) {
        tf = T_flip_pro_idx(pts, triads, slump, ids);
        nits = (int)ids.size();
        nit++;
        if (tf < 0)
        {
            throw Exception("cannot triangualte this set ");

            return -4;
        }
    }
    
    return 1;
}

void circle_cent4(TBase r1, TBase c1, TBase r2, TBase c2, TBase r3, TBase c3, TBase& r, TBase& c, TBase& ro2) {
    /*
    *  function to return the center of a circle and its radius
    * degenerate case should never be passed to this routine!!!!!!!!!!!!!
    * but will return r0 = -1 if it is.
    */

    TBase rd, cd;
    TBase v1 = 2 * (r2 - r1);
    TBase v2 = 2 * (c2 - c1);
    TBase v3 = r2*r2 - r1*r1 + c2*c2 - c1*c1;
    TBase v4 = 2 * (r3 - r1);
    TBase v5 = 2 * (c3 - c1);
    TBase v6 = r3*r3 - r1*r1 + c3*c3 - c1*c1;
    TBase v7 = v2*v4 - v1*v5;
    
    if (v7 == 0) {
        r = 0;
        c = 0;
        ro2 = -1;
        return;
    }

    cd = (v4*v3 - v1*v6) / v7;
    if (v1 != 0)
        rd = (v3 - c*v2) / v1;
    else
        rd = (v6 - c*v5) / v4;

    ro2 = (TBase)((rd - r1)*(rd - r1) + (cd - c1)*(cd - c1));
    r = (TBase)rd;
    c = (TBase)cd;
}

/*
flip pairs of triangles that are not valid delaunay triangles
the Cline-Renka test is used rather than the less stable circum
circle center computation test of s-hull.

or the more expensive determinant test.

*/
int T_flip_pro(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, int numt, int start, std::set<int>& ids) {
    TBase r3, c3;
    int pa, pb, pc, pd, D, L1, L2, L3, L4, T2;

    Triad tx, tx2;

    for (int t = start; t < numt; ++t) {
        Triad &tri = triads[t];
        // test all 3 neighbours of tri 

        int flipped = 0;

        if (tri.bc >= 0) {
            pa = slump[tri.a];
            pb = slump[tri.b];
            pc = slump[tri.c];

            T2 = tri.bc;
            Triad &t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t){
                D = t2.c;
                pd = slump[t2.c];

                if (tri.b == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.b == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.b == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            if (pd < 0 || pd > 100)
                int dfx = 9;

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pa].r, pts[pa].c, pts[pb].r, pts[pb].c, pts[pc].r, pts[pc].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;
                L2 = tri.ac;
                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.a;
                    tx.b = tri.b;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;


                    // triangle 2;
                    tx2.a = tri.a;
                    tx2.b = tri.c;
                    tx2.c = D;

                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;


                    ids.insert(t);
                    ids.insert(T2);

                    t2 = tx2;
                    tri = tx;
                    flipped = 1;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad &t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad &t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }
                }
            }
        }

        if (flipped == 0 && tri.ab >= 0){
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ab;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pc].r, pts[pc].c, pts[pb].r, pts[pb].c, pts[pa].r, pts[pa].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ac;
                L2 = tri.bc;
                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.c;
                    tx.b = tri.a;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;


                    // triangle 2;
                    tx2.a = tri.c;
                    tx2.b = tri.b;
                    tx2.c = D;

                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;


                    ids.insert(t);
                    ids.insert(T2);

                    t2 = tx2;
                    tri = tx;
                    flipped = 1;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad& t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad& t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }
                }
            }
        }


        if (flipped == 0 && tri.ac >= 0) {
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ac;
            Triad &t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t){
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a){
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t){
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b){
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pb].r, pts[pb].c, pts[pa].r, pts[pa].c, pts[pc].r, pts[pc].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;   // .ac shared limb
                L2 = tri.bc;
                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.b;
                    tx.b = tri.a;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;


                    // triangle 2;
                    tx2.a = tri.b;
                    tx2.b = tri.c;
                    tx2.c = D;

                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;

                    ids.insert(t);
                    ids.insert(T2);

                    t2 = tx2;
                    tri = tx;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad& t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad& t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }
                }
            }
        }
    }

    return 1;
}

int Cline_Renka_test(TBase& Ax, TBase& Ay,
    TBase& Bx, TBase& By,
    TBase& Cx, TBase& Cy,
    TBase& Dx, TBase& Dy)
{
    TBase v1x = Bx - Ax;
    TBase v1y = By - Ay;
    TBase v2x = Cx - Ax;
    TBase v2y = Cy - Ay;
    TBase v3x = Bx - Dx;
    TBase v3y = By - Dy;
    TBase v4x = Cx - Dx;
    TBase v4y = Cy - Dy;
    TBase cosA = v1x*v2x + v1y*v2y;
    TBase cosD = v3x*v4x + v3y*v4y;

    if (cosA < 0 && cosD < 0) // two obtuse angles 
        return -1 ;

    TBase ADX = Ax - Dx;
    TBase ADy = Ay - Dy;

    if (cosA > 0 && cosD > 0)  // two acute angles
        return 1;

    TBase sinA = fabs(v1x*v2y - v1y*v2x);
    TBase sinD = fabs(v3x*v4y - v3y*v4x);

    if (cosA*sinD + sinA*cosD < 0)
        return -1;

    return 1;
}

// same again but with set of triangle ids to be iterated over.
int T_flip_pro_idx(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, std::set<int>& ids) {
    TBase r3, c3;
    int pa, pb, pc, pd, D, L1, L2, L3, L4, T2;

    Triad tx, tx2;
    std::set<int> ids2;
    ids2.clear();

    std::set<int>::const_iterator x = ids.begin();
    while (x != ids.end()){
        int t = *x;
        x++;

        Triad& tri = triads[t];
        // test all 3 neighbours of tri 
        int flipped = 0;

        if (tri.bc >= 0) {
            pa = slump[tri.a];
            pb = slump[tri.b];
            pc = slump[tri.c];

            T2 = tri.bc;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.b == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.b == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.b == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error ");
                return -6;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pa].r, pts[pa].c, pts[pb].r, pts[pb].c, pts[pc].r, pts[pc].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;
                L2 = tri.ac;

                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.a;
                    tx.b = tri.b;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;

                    // triangle 2;
                    tx2.a = tri.a;
                    tx2.b = tri.c;
                    tx2.c = D;

                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;

                    ids2.insert(t);
                    ids2.insert(T2);

                    t2 = tx2;
                    tri = tx;
                    flipped = 1;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad& t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad& t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }
                }
            }
        }

        if (flipped == 0 && tri.ab >= 0) {
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ab;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -6;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pc].r, pts[pc].c, pts[pb].r, pts[pb].c, pts[pa].r, pts[pa].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ac;
                L2 = tri.bc;
                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.c;
                    tx.b = tri.a;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;


                    // triangle 2;
                    tx2.a = tri.c;
                    tx2.b = tri.b;
                    tx2.c = D;

                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;


                    ids2.insert(t);
                    ids2.insert(T2);

                    t2 = tx2;
                    tri = tx;
                    flipped = 1;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad& t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad& t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }

                }
            }
        }

        if (flipped == 0 && tri.ac >= 0) {
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ac;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            }
            else{
                throw Exception("triangle flipping error. ");
                return -6;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pb].r, pts[pb].c, pts[pc].r, pts[pc].c, pts[pa].r, pts[pa].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;   // .ac shared limb
                L2 = tri.bc;
                if (L1 != L3 && L2 != L4) {  // need this check for stability.
                    tx.a = tri.b;
                    tx.b = tri.a;
                    tx.c = D;

                    tx.ab = L1;
                    tx.ac = T2;
                    tx.bc = L3;


                    // triangle 2;
                    tx2.a = tri.b;
                    tx2.b = tri.c;
                    tx2.c = D;


                    tx2.ab = L2;
                    tx2.ac = t;
                    tx2.bc = L4;


                    ids2.insert(t);
                    ids2.insert(T2);

                    t2 = tx2;
                    tri = tx;

                    // change knock on triangle labels.
                    if (L3 >= 0) {
                        Triad& t3 = triads[L3];
                        if (t3.ab == T2)
                            t3.ab = t;
                        else if (t3.bc == T2)
                            t3.bc = t;
                        else if (t3.ac == T2)
                            t3.ac = t;
                    }

                    if (L2 >= 0) {
                        Triad& t4 = triads[L2];
                        if (t4.ab == t)
                            t4.ab = T2;
                        else if (t4.bc == t)
                            t4.bc = T2;
                        else if (t4.ac == t)
                            t4.ac = T2;
                    }
                }
            }
        }
    }

    ids.clear();
    ids.insert(ids2.begin(), ids2.end());

    return 1;
}

/* test the seed configuration to see if the center
of the circum circle lies inside the seed triangle.

if not issue a warning.
*/
int test_center(Shx& pt0, Shx& pt1, Shx& pt2) {
    TBase r01 = pt1.r - pt0.r;
    TBase c01 = pt1.c - pt0.c;

    TBase r02 = pt2.r - pt0.r;
    TBase c02 = pt2.c - pt0.c;

    TBase r21 = pt1.r - pt2.r;
    TBase c21 = pt1.c - pt2.c;

    TBase v = r01*r02 + c01*c02;
    if (v < 0)
        return -1;

    v = r21*r02 + c21*c02;
    if (v > 0)
        return -1;

    v = r01*r21 + c01*c21;
    if (v < 0)
        return -1;

    return 1;
}

int T_flip_edge(std::vector<Shx>& pts, std::vector<Triad>& triads, std::vector<int>& slump, int numt, int start, std::set<int>& ids) {
    TBase r3, c3;
    int pa, pb, pc, pd, D, L1, L2, L3, L4, T2;

    Triad tx, tx2;
    for (int t = start; t < numt; ++t){
        Triad &tri = triads[t];
        // test all 3 neighbours of tri 

        int flipped = 0;
        if (tri.bc >= 0 && (tri.ac < 0 || tri.ab < 0)) {
            pa = slump[tri.a];
            pb = slump[tri.b];
            pc = slump[tri.c];

            T2 = tri.bc;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.b == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            } else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.b == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.b == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            if (pd < 0 || pd > 100)
                int dfx = 9;

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pa].r, pts[pa].c, pts[pb].r, pts[pb].c, pts[pc].r, pts[pc].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;
                L2 = tri.ac;
                //	if( L1 != L3 && L2 != L4 ){  // need this check for stability.

                tx.a = tri.a;
                tx.b = tri.b;
                tx.c = D;

                tx.ab = L1;
                tx.ac = T2;
                tx.bc = L3;


                // triangle 2;
                tx2.a = tri.a;
                tx2.b = tri.c;
                tx2.c = D;

                tx2.ab = L2;
                tx2.ac = t;
                tx2.bc = L4;


                ids.insert(t);
                ids.insert(T2);

                t2 = tx2;
                tri = tx;
                flipped = 1;

                // change knock on triangle labels.
                if (L3 >= 0) {
                    Triad& t3 = triads[L3];
                    if (t3.ab == T2)
                        t3.ab = t;
                    else if (t3.bc == T2)
                        t3.bc = t;
                    else if (t3.ac == T2)
                        t3.ac = t;
                }

                if (L2 >= 0) {
                    Triad& t4 = triads[L2];
                    if (t4.ab == t)
                        t4.ab = T2;
                    else if (t4.bc == t)
                        t4.bc = T2;
                    else if (t4.ac == t)
                        t4.ac = T2;
                }
                //	}
            }
        }


        if (flipped == 0 && tri.ab >= 0 && (tri.ac < 0 || tri.bc < 0)) {
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ab;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t){
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pc].r, pts[pc].c, pts[pb].r, pts[pb].c, pts[pa].r, pts[pa].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ac;
                L2 = tri.bc;
                //	if( L1 != L3 && L2 != L4 ){  // need this check for stability.

                tx.a = tri.c;
                tx.b = tri.a;
                tx.c = D;

                tx.ab = L1;
                tx.ac = T2;
                tx.bc = L3;


                // triangle 2;
                tx2.a = tri.c;
                tx2.b = tri.b;
                tx2.c = D;

                tx2.ab = L2;
                tx2.ac = t;
                tx2.bc = L4;


                ids.insert(t);
                ids.insert(T2);

                t2 = tx2;
                tri = tx;
                flipped = 1;

                // change knock on triangle labels.
                if (L3 >= 0) {
                    Triad& t3 = triads[L3];
                    if (t3.ab == T2)
                        t3.ab = t;
                    else if (t3.bc == T2)
                        t3.bc = t;
                    else if (t3.ac == T2)
                        t3.ac = t;
                }

                if (L2 >= 0) {
                    Triad& t4 = triads[L2];
                    if (t4.ab == t)
                        t4.ab = T2;
                    else if (t4.bc == t)
                        t4.bc = T2;
                    else if (t4.ac == t)
                        t4.ac = T2;
                }

                //	}

            }
        }


        if (flipped == 0 && tri.ac >= 0 && (tri.bc < 0 || tri.ab < 0)) {
            pc = slump[tri.c];
            pb = slump[tri.b];
            pa = slump[tri.a];

            T2 = tri.ac;
            Triad& t2 = triads[T2];
            // find relative orientation (shared limb).
            if (t2.ab == t) {
                D = t2.c;
                pd = slump[t2.c];

                if (tri.a == t2.a) {
                    L3 = t2.ac;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ac;
                }
            }
            else if (t2.ac == t) {
                D = t2.b;
                pd = slump[t2.b];

                if (tri.a == t2.a) {
                    L3 = t2.ab;
                    L4 = t2.bc;
                } else {
                    L3 = t2.bc;
                    L4 = t2.ab;
                }
            }
            else if (t2.bc == t) {
                D = t2.a;
                pd = slump[t2.a];

                if (tri.a == t2.b) {
                    L3 = t2.ab;
                    L4 = t2.ac;
                } else {
                    L3 = t2.ac;
                    L4 = t2.ab;
                }
            } else {
                throw Exception("triangle flipping error. ");
                return -5;
            }

            r3 = pts[pd].r;
            c3 = pts[pd].c;

            int XX = Cline_Renka_test(pts[pb].r, pts[pb].c, pts[pa].r, pts[pa].c, pts[pc].r, pts[pc].c, r3, c3);

            if (XX < 0) {
                L1 = tri.ab;   // .ac shared limb
                L2 = tri.bc;
                //	if( L1 != L3 && L2 != L4 ){  // need this check for stability.

                tx.a = tri.b;
                tx.b = tri.a;
                tx.c = D;

                tx.ab = L1;
                tx.ac = T2;
                tx.bc = L3;


                // triangle 2;
                tx2.a = tri.b;
                tx2.b = tri.c;
                tx2.c = D;

                tx2.ab = L2;
                tx2.ac = t;
                tx2.bc = L4;

                ids.insert(t);
                ids.insert(T2);

                t2 = tx2;
                tri = tx;

                // change knock on triangle labels.
                if (L3 >= 0) {
                    Triad& t3 = triads[L3];
                    if (t3.ab == T2)
                        t3.ab = t;
                    else if (t3.bc == T2)
                        t3.bc = t;
                    else if (t3.ac == T2)
                        t3.ac = t;
                }

                if (L2 >= 0) {
                    Triad& t4 = triads[L2];
                    if (t4.ab == t)
                        t4.ab = T2;
                    else if (t4.bc == t)
                        t4.bc = T2;
                    else if (t4.ac == t)
                        t4.ac = T2;
                }

                //}
            }
        }
    }

    return 1;
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

    vector<Integers> graph;
    graph.resize(m, vector<int>());
    vector<bool> used;
    used.resize(m);
    queue<int> qu;

    std::vector<Shx> pts(m);
    for (int i = 0; i < m; ++i) {
        pts[i].r = dpoints[i]._x;
        pts[i].c = dpoints[i]._y;
        pts[i].id = i;
    }

    std::vector<Triad> triads;
    if (m > 2) {
        s_hull_pro(pts, triads);
    }

    for (size_t i = 0; i < triads.size(); ++i) {
        Triad& t = triads[i];
        graph[t.a].push_back(t.b);
        graph[t.b].push_back(t.a);
        graph[t.a].push_back(t.c);
        graph[t.c].push_back(t.a);
        graph[t.c].push_back(t.b);
        graph[t.b].push_back(t.c);
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
            const Triad& t = triads[0];
            temp.push_back(t.a);
            temp.push_back(t.b);
            temp.push_back(t.c);
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