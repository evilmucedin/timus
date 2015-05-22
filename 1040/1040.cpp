#define _CRT_SECURE_NO_WARNINGS

#include    <cstdio>
#include	<cstdlib>

#include	<vector>
#include	<queue>

using namespace std;

typedef vector<int> Integers;
typedef vector<Integers> Graph;

int main(void)
{
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif
	int n, m;
	scanf("%d %d", &n, &m);

	Graph g(n);
	Integers start(m);
	Integers finish(m);
	for (int i = 0; i < m; i++)
	{
		int a, b;
		scanf("%d %d", &a, &b);
		a--;
		b--;
		g[a].push_back(i);
		g[b].push_back(i);
		start[i] = a;
		finish[i] = b;
	}

	Graph g2(m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < g[i].size(); ++j) {
			for (int k = j + 1; k < g[i].size(); ++k) {
				g2[g[i][j]].push_back( g[i][k] );
				g2[g[i][k]].push_back(g[i][j]);
			}
		}
	}

	Integers fn(m, -1);
	queue<int> q;
	int index = 1;
	q.push(0);
	while (!q.empty()) {
		int v = q.front();
		q.pop();
		for (int i = 0; i < g[v].size(); ++i) {
			if (fn[g[v][i]] == -1) {
				fn[g[v][i]] = index++;
				q.push(start[g[v][i]]);
				q.push(finish[g[v][i]]);
			}
		}
	}

	printf("YES\n");
	for (int i = 0; i < m; i++)
		printf("%d ", fn[i]);
	printf("\n");

	return 0;
}

