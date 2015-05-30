#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>

using namespace std;

struct TEdge {
	int _begin;
	int _end;
};

typedef vector<TEdge> TEdges;
typedef vector<int> TIntegers;

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif
	
	int n;
	int m;
	scanf("%d%d", &n, &m);

	TEdges edges(m);
	for (int i = 0; i < m; ++i) {
		scanf("%d%d", &edges[i]._begin, &edges[i]._end);
		--edges[i]._begin;
		--edges[i]._end;
	}

	TIntegers order(n);
	for (int i = 0; i < n; ++i) {
		scanf("%d", &order[i]);
		--order[i];
	}

	TIntegers pos(n);
	for (int i = 0; i < n; ++i) {
		pos[order[i]] = i;
	}

	int i = 0;
	while (i < edges.size() && pos[edges[i]._begin] < pos[edges[i]._end])
		++i;

	printf((i == edges.size()) ? "YES\n" : "NO\n");

	return 0;
}