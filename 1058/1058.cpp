#include	<cstdio>
#include	<cmath>

#define M_PI 3.14159265358979323846

#define	MAXP	200
#define	INF		1E10
#define	EPS		1E-9
#define	DA		1E-6

#define	SQR(X)  ((X)*(X))

struct TPoint
{
	double x;
	double y;
};

int n;
TPoint points[MAXP];
TPoint center;
double minlenght = INF;
double bestangle;
double S0;

inline double getDist(TPoint p1, TPoint p2)
{
	return sqrt(SQR(p1.x - p2.x) + SQR(p1.y - p2.y));
}

inline double getTriangleSquare(TPoint p1, TPoint p2, TPoint p3)
{
	return fabs((p2.x - p3.x)*(p1.y - p3.y) - (p2.y - p3.y)*(p1.x - p3.x)) / 2.0;
}

inline bool intersect(double angle, TPoint p1, TPoint p2, TPoint* p)
{
	double A1 = cos(angle);
	double B1 = p1.x - p2.x;
	double C1 = p1.x - center.x;
	double A2 = sin(angle);
	double B2 = p1.y - p2.y;
	double C2 = p1.y - center.y;
	double D = A1*B2 - A2*B1;
	if (fabs(D) > EPS)
	{
		double u = (C1*B2 - B1*C2) / D;
		double v = (A1*C2 - C1*A2) / D;
		p->x = p1.x + (p2.x - p1.x)*v;
		p->y = p1.y + (p2.y - p1.y)*v;
		return (u > -EPS) && (v > -EPS) && (v < 1.0 + EPS);
	}
	else
	{
		return false;
	}
}

inline void getPoint(double angle, TPoint* point, int* side)
{
	while (!intersect(angle, points[*side], points[((*side) + 1) % n], point))
		*side = ((*side) + 1) % n;
}

int side = 0;

inline double Length(double angle)
{
	TPoint point;
	getPoint(angle, &point, &side);

	double S = 0.0;
	double TS;
	int index = side;
	while (2.0*S < S0)
	{
		TS = getTriangleSquare(point, points[(index + 1) % n], points[(index + 2) % n]);
		index = (index + 1) % n;
		S += TS;
	}

	double k = 1.0 - (S - S0 / 2.0) / TS;

	TPoint point2;
	point2.x = points[index].x + (points[(index + 1) % n].x - points[index].x)*k;
	point2.y = points[index].y + (points[(index + 1) % n].y - points[index].y)*k;

	//double SS = getTriangleSquare(points[2], point, point2);

	return getDist(point, point2);
}

inline double TryLength(double angle)
{
	double length = Length(angle);
	if (length < minlenght)
	{
		minlenght = length;
		bestangle = angle;
	}
	return length;
}

inline double getAngle(double dy, double dx)
{
	return atan2(dy, dx);
	/*
	if (fabs(dx) > EPS)
	{
		double angle = atan(fabs(dy / dx));
		if ((dy > 0) && (dx > 0))
			return angle;
		if ((dy > 0) && (dx <= 0))
			return M_PI - angle;
		if ((dy <= 0) && (dx <= 0))
			return M_PI + angle;
		return 2.0*M_PI - angle;
	}
	else
	{
		if (dy > 0.0)
			return M_PI / 2.0;
		else
			return 3.0*M_PI / 2.0;
	}
	*/
}

inline double FixAngle(double angle)
{
	double len1 = TryLength(angle - DA);
	double len2 = TryLength(angle);
	double len3 = TryLength(angle + DA);
	double fs = (len3 - len1) / 2.0 / DA;
	double fss = (len3 - 2.0*len2 + len1) / DA / DA;
	if (fabs(fss) > 1E-15)
		angle -= fs / fss;
	return angle;
}

void Solve(void)
{
	double angle;
	for (angle = 0.0; angle < 2.0*M_PI + EPS; angle += 0.0001)
	{
		TryLength(angle);
		double angle1 = angle;
		for (int i = 0; i < 3; ++i)
			angle1 = FixAngle(angle1);
	}

	int i;
	for (i = 0; i < n; ++i)
	{
		double angle1 = getAngle(points[i].y - center.y, points[i].x - center.x);
		double angle2 = getAngle(points[(i + 1) % n].y - center.y, points[(i + 1) % n].x - center.x);
		static const int N = 1000;
		for (int j = 0; j < N; j++)
		{
			double angle = (angle1*j + angle2*(N - j)) / N;
			TryLength(angle);
			for (int i = 0; i < 3; ++i)
				angle = FixAngle(angle);
		}
	}

	for (int j = 0; j < 10; ++j)
	{
		double angle = bestangle;
		for (int i = 0; i < 10; ++i)
			angle = FixAngle(angle);
	}
}

int main(void)
{
	scanf("%d", &n);
	int i;
	for (i = 0; i < n; ++i)
		scanf("%lf %lf", &(points[i].x), &(points[i].y));

	double SX = 0.0;
	double SY = 0.0;
	for (i = 0; i < n; ++i)
		SX += points[i].x;
	for (i = 0; i < n; ++i)
		SY += points[i].y;
	center.x = SX / n;
	center.y = SY / n;

	S0 = 0.0;
	for (i = 1; i < n - 1; ++i)
		S0 += getTriangleSquare(points[0], points[i], points[i + 1]);

	if (fabs(S0) > EPS)
	{
		Solve();
		printf("%.4lf\n", minlenght);
	}
	else
	{
		printf("%.4lf\n", 0.0);
	}

	return 0;
}

