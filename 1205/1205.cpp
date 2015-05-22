#include	<cstdio>
#include    <cmath>

#include	<iostream>
#include	<fstream>

using namespace std;

#define	MAXN	205

struct CPoint
{
	double x, y;
};

typedef double CLine[MAXN];
typedef double CGraph[MAXN][MAXN];
typedef	CPoint CPoints[MAXN];
typedef int CILine[MAXN];
typedef int CIGraph[MAXN][MAXN];

int n;
double velfoot, velmetro;

CGraph graph;
CLine time;
CPoints coords;
CILine previous;
CILine F;
CIGraph conn;

void Input(void)
{
	cin >> velfoot >> velmetro;
	cin >> n;

	for (int i = 0; i < n; i++)
		cin >> coords[i].x >> coords[i].y;

	for (int i = 0; i < n + 2; i++)
		for (int j = 0; j < n + 2; j++)
			conn[i][j] = 0;

	int a, b;
	cin >> a >> b;
	while ((a != 0) || (b != 0))
	{
		conn[a - 1][b - 1] = 1;
		conn[b - 1][a - 1] = 1;
		cin >> a >> b;
	}

	cin >> coords[n].x >> coords[n].y;
	cin >> coords[n + 1].x >> coords[n + 1].y;
}

void MakeGraph(void)
{
	int i, j;
	for (i = 0; i < n + 2; i++)
		for (j = 0; j < n + 2; j++)
			graph[i][j] = sqrt((coords[i].x - coords[j].x) * (coords[i].x - coords[j].x) + (coords[i].y - coords[j].y) * (coords[i].y - coords[j].y)) / velfoot;
	for (i = 0; i < n + 2; i++)
		for (j = 0; j < n + 2; j++)
			if (conn[i][j])
				graph[i][j] = sqrt((coords[i].x - coords[j].x) * (coords[i].x - coords[j].x) + (coords[i].y - coords[j].y) * (coords[i].y - coords[j].y)) / velmetro;

	/*for (i = 0; i < n + 2; i++)
	{
	for (j = 0; j < n + 2; j++)
	cout << graph[i][j] << " ";
	cout << "\n";
	}*/
}

void Dijkstra(void)
{
	time[n] = 0;
	int i;
	for (i = 0; i < n + 2; i++)
		if (i != n)
		{
			time[i] = graph[n][i];
			previous[i] = n;
			F[i] = 1;
		}
		else
			F[i] = 0;

	for (i = 0; i < n + 1; i++)
	{
		int w, k;
		double min = 1E10;
		for (k = 0; k < n + 2; k++)
			if (F[k] && (time[k] < min))
			{
				min = time[k];
				w = k;
			}
		F[w] = 0;

		for (k = 0; k < n + 2; k++)
			if (F[k])
				if (time[w] + graph[w][k] < time[k])
				{
					time[k] = time[w] + graph[w][k];
					previous[k] = w;
				}
	}
}

void Output(void)
{
	//for (int l = 0; l < n + 2; l++)
	//	cout << time[l] << " ";
	//cout << "\n";

	CILine path;
	int pathc = 0;
	int now = n + 1;
	while (now != n)
	{
		path[pathc] = now;
		now = previous[now];
		pathc++;
	}

	printf("%.7f\n", time[n + 1]);
	printf("%d ", pathc - 1);
	for (int i = pathc - 1; i >= 1; i--)
		printf("%d ", path[i] + 1);
}

int main(void)
{
	Input();
	MakeGraph();
	Dijkstra();
	Output();
	return 0;
}
