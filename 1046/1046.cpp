#define _CRT_SECURE_NO_WARNINGS

#include <cmath>

#include <iostream>
#include <vector>

using namespace std;

struct TPoint {
	long double _x;
	long double _y;

	TPoint() {

	}

	TPoint(long double x, long double y)
		: _x(x)
		, _y(y)
	{

	}
};
typedef vector<TPoint> TPoints;
typedef vector<long double> TLongDoubles;

void Rotate(const TPoint& center, const TPoint& p, long double alpha, TPoint* result) {
	long double dx = p._x - center._x;
	long double dy = p._y - center._y;
	result->_x = center._x + dx*cos(alpha) - dy*sin(alpha);
	result->_y = center._y + dx*sin(alpha) + dy*cos(alpha);
}

void RotateAll(const TPoint& p, const TPoints& points, const TLongDoubles& angles, TPoint* result) {
	*result = p;
	for (size_t i = 0; i < points.size(); ++i) {
		Rotate(points[i], *result, angles[i], result);
	}
}

void Out(long double x) {
	if (x < 0) {
		printf("-");
		x = -x;
	}
	int ix = x*100.0;
	printf("%d", ix/100);
	if (ix % 100) {
		printf(".%d", (ix % 100)/10);
		if ((ix % 100) % 10) {
			printf("%d", (ix % 100) % 10);
		}
	}
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	cin >> n;

	TPoints points(n);
	for (int i = 0; i < n; ++i) {
		scanf("%Lf%Lf", &points[i]._x, &points[i]._y);
	}

	TLongDoubles angles(n);
	for (int i = 0; i < n; ++i) {
		scanf("%Lf", &angles[i]);
		angles[i] *= 2.0*acos(0.0) / 180.0;
	}

	TPoint p1(0.0, 1.0);
	TPoint res1;
	RotateAll(p1, points, angles, &res1);
	TPoint p2(1.0, 0.0);
	TPoint res2;
	RotateAll(p2, points, angles, &res2);

	TPoint cn1((p1._x + res1._x)/2.0, (p1._y + res1._y) / 2.0);
	TPoint cn2((p2._x + res2._x) / 2.0, (p2._y + res2._y) / 2.0);

	long double a1 = p1._x - res1._x;
	long double b1 = p1._y - res1._y;
	long double c1 = -a1*cn1._x -b1*cn1._y;
	long double a2 = p2._x - res2._x;
	long double b2 = p2._y - res2._y;
	long double c2 = -a2*cn2._x - b2*cn2._y;

	long double d = a1*b2 - a2*b1;
	long double u = (c1*b2 - b1*c2) / d;
	long double v = (a1*c2 - c1*a2) / d;
	TPoint c(-u, -v);

	for (int i = 0; i < n; ++i) {
		Out(c._x);
		printf(" ");
		Out(c._y);
		printf("\n");
		Rotate(points[i], c, angles[i], &c);
	}

	return 0;
}