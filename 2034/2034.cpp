#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>
#include <queue>

using namespace std;

typedef vector<int> TLine;
typedef vector<TLine> TGraph;
typedef queue<int> TQueue;

void DFS(const TGraph& g, int v, TLine* result) {
	result->resize(g.size());
	fill(result->begin(), result->end(), -1);
	TQueue q;
	(*result)[v] = 0;
	q.push(v);
	while (!q.empty()) {
		int now = q.front();
		q.pop();
		for (int i = 0; i < g[now].size(); ++i) {
			int next = g[now][i];
			if (-1 == (*result)[next]) {
				(*result)[next] = (*result)[now] + 1;
				q.push(next);
			}
		}
	}
}

void TrackMaxDist(const TGraph& g, const TLine& sDist, const TLine& rDist, int now, TLine* cache) {
	if (-1 == (*cache)[now]) {
		int prevs = -1;
		for (int i = 0; i < g[now].size(); ++i) {
			int next = g[now][i];
			if (sDist[next] + 1 == sDist[now]) {
				TrackMaxDist(g, sDist, rDist, next, cache);
                prevs = max(prevs, (*cache)[next]);
			}
		}
        int mx;
        if (-1 != prevs) {
            mx = min(prevs, rDist[now]);
        } else {
            mx = rDist[now];
        }
		(*cache)[now] = mx;		
	}
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	int m;
	scanf("%d%d", &n, &m);

	TGraph g;
	g.resize(n);

	for (int i = 0; i < m; ++i) {
		int begin;
		int end;
		scanf("%d%d", &begin, &end);
		--begin;
		--end;
		g[begin].push_back(end);
		g[end].push_back(begin);
	}

	int s;
	int f;
	int r;
	scanf("%d%d%d", &s, &f, &r);
	--s;
	--f;
	--r;

	TLine sDist;
	DFS(g, s, &sDist);
	TLine rDist;
	DFS(g, r, &rDist);

	TLine cache(n, -1);
	cache[s] = rDist[s];

	TrackMaxDist(g, sDist, rDist, f, &cache);

	printf("%d\n", cache[f]);

	return 0;
}