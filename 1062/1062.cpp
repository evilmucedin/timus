#define _CRT_SECURE_NO_WARNINGS

#include	<cstdio>
#include	<cmath>
#include	<vector>

using namespace std;

#define	MAXN	150
#define	EPS		1E-15

class TPoint
{
public:
	double x;
	double y;
	double z;
};

int n;
TPoint v[MAXN];

void Input(void)
{
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%lf %lf %lf", &(v[i].x), &(v[i].y), &(v[i].z));
}

typedef vector<TPoint> TPoints;

inline double Comb(TPoint& a, TPoint& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

#define	SQR(X)  ((X)*(X))

inline double Dist(TPoint& a, TPoint& b)
{
	return sqrt(SQR(a.x - b.x) + SQR(a.y - b.y) + SQR(a.z - b.z));
}

inline void Intersect(TPoint& a, TPoint& b, TPoint& c, TPoint& res)
{
	double t = -(c.x*a.x + c.y*a.y + c.z*a.z) / (c.x*(b.x - a.x) + c.y*(b.y - a.y) + c.z*(b.z - a.z));
	res.x = a.x + t*(b.x - a.x);
	res.y = a.y + t*(b.y - a.y);
	res.z = a.z + t*(b.z - a.z);
}

int NotEmpty(int p)
{
	TPoints current;
	TPoint temp;
	temp.x = 1.0;
	temp.y = 0.0;
	temp.z = 0.0;
	current.push_back(temp);
	temp.x = 0.0;
	temp.y = 1.0;
	temp.z = 0.0;
	current.push_back(temp);
	temp.x = 0.0;
	temp.y = 0.0;
	temp.z = 1.0;
	current.push_back(temp);

	int i;
	for (i = 0; i < n; i++)
		if (i != p)
		{
			TPoints nw;
			TPoint c;
			int j;
			c.x = (v[p].x - v[i].x) / (v[p].x*v[i].x);
			c.y = (v[p].y - v[i].y) / (v[p].y*v[i].y);
			c.z = (v[p].z - v[i].z) / (v[p].z*v[i].z);

			j = 0;
			while ((j < current.size()) && (Comb(c, current[j]) <= EPS))
				j++;
			if (j == current.size())
				break;

			int r = 0;
			while (r < current.size())
			{
				if (Comb(c, current[j]) > EPS)
				{
					nw.push_back(current[j]);
					j = (j + 1) % current.size();
					r++;
				}
				else
				{
					int k = (j + 1) % current.size();
					int l = 0;
					r++;
					while ((l < current.size() + 1) && (Comb(c, current[k]) <= EPS))
					{
						k = (k + 1) % current.size();
						l++;
						r++;
					}
					if (l < current.size() + 1)
					{
						TPoint nwPoint;
						Intersect(nw[nw.size() - 1], current[j], c, nwPoint);
						nw.push_back(nwPoint);
						Intersect(current[(k + current.size() - 1) % current.size()], current[k], c, nwPoint);
						nw.push_back(nwPoint);
					}
					j = k;
				}
			}

			/*for (int i = 0; i < nw.size(); i++)
			printf("    %.3lf %.3lf %.3lf\n", nw[i].x, nw[i].y, nw[i].z);
			printf("\n\n");*/

			current = nw;
		}

	if (i != n)
	{
		return 0;
	}
	else
	{
		double d = -1.0;
		for (int i = 0; i < current.size() - 1; i++)
			for (int j = i + 1; j < current.size() - 1; j++)
				if (Dist(current[i], current[j]) > d)
					d = Dist(current[i], current[j]);
		return d > EPS;
	}
}

int main(void)
{
	/*TPoint a, b, c, res;
	a.x = 1.0;
	a.y = 0.0;
	a.z = 0.0;
	b.x = -1.0;
	b.y = 1.0;
	b.z = 0.0;
	c.x = 1.0;
	c.y = 0.0;
	c.z = 0.0;
	Intersect(a, b, c, res);*/
	Input();
	for (int i = 0; i < n; i++)
	{
		if (NotEmpty(i))
			printf("Yes\n");
		else
			printf("No\n");
	}

	return 0;
}
