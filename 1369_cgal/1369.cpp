#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

using namespace std;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation_2;

typedef vector<Point_2> TPoints;

int main()
{	
	int m;	
	cin >> m;
	
	TPoints points(m);

	for (int i = 0; i < m; ++i)
	{
		long double x;
		long double y;
		cin >> x >> y;
		points.push_back(Point_2(x, y));
	}


	Delaunay_triangulation_2 dt(points.begin(), points.end());
	
	int n;
	cin >> n;
	for (int i = 0; i < n; ++i)
	{
		long double x;
		long double y;
		cin >> x >> y;
		Point_2 p(x, y);
		dt.nearest_vertex(p);
	}


	return 0;
}
