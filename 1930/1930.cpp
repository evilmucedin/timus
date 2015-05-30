#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>
#include <queue>
using namespace std;

struct TEdge {
	int _end;
	int _cost;
	
	TEdge() {

	}

	TEdge(int end, int cost)
		: _end(end)
		, _cost(cost)
	{

	}
};

typedef vector<TEdge> TEdges;
typedef vector<TEdges> TGraph;

struct TQueueItem {
	int _v;
	int _cost;

	TQueueItem() {

	}

	TQueueItem(int v, int cost) 
		: _v(v)
		, _cost(cost)
	{
	
	}

	bool operator<(const TQueueItem& q) const {
		return _cost > q._cost;
	}
};

typedef priority_queue<TQueueItem> TQueue;

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	int m;
	scanf("%d%d", &n, &m);

	int nn = 2 * n;
	TGraph g(nn);
	for (int i = 0; i < n; ++i) {
		g[2 * i].push_back( TEdge(2*i + 1, 1) );
		g[2 * i + 1].push_back(TEdge(2 * i, 1));
	}

	for (int i = 0; i < m; ++i) {
		int begin;
		int end;
		scanf("%d%d", &begin, &end);
		--begin;
		--end;
		g[2 * begin].push_back( TEdge(2*end, 0) );
		g[2 * end + 1].push_back(TEdge(2 * begin + 1, 0));
	}

	int start;
	int finish;
	scanf("%d%d", &start, &finish);
	--start;
	--finish;
	g[2 * start].push_back(TEdge(2 * start + 1, 0));
	g[2 * start + 1].push_back(TEdge(2 * start, 0));

	static const int INF = 123456789;

	TQueue q;
	vector<int> dist(nn);
	for (int i = 0; i < nn; ++i) {
		dist[i] = INF;
	}
	vector<bool> visited(nn);
	int iStart = 2 * start;
	dist[iStart] = 0;
	q.push(TQueueItem(iStart, 0));
	while (!q.empty()) {
		TQueueItem item = q.top();
		q.pop();
		if (!visited[item._v]) {
			visited[item._v] = true;
			for (TEdges::const_iterator toEdge = g[item._v].begin(); toEdge != g[item._v].end(); ++toEdge) {
				int newDist = dist[item._v] + toEdge->_cost;
				if (dist[toEdge->_end] > newDist) {
					dist[toEdge->_end] = newDist;
					q.push(TQueueItem(toEdge->_end, newDist));
				}
			}
		}
	}
	int ans = dist[2 * finish];
	ans = min(ans, dist[2 * finish + 1]);

	printf("%d\n", ans);

	return 0;
}