#include <cstdio>
#include <cmath>

#include <vector>
#include <queue>

using namespace std;

struct TNode {
	int _to;
	double _value;
	double _p;

	TNode(int to, double value, double p)
		: _to(to)
		, _value(value)
		, _p(p)
	{

	}

	TNode() {
	}
};

struct TCost {
	int _iDist;
	double _fine;

	TCost(int iDist, double fine)
		: _iDist(iDist)
		, _fine(fine)
	{

	}

	TCost()
	{

	}
};

typedef vector<TNode> TNodes;
typedef vector<TNodes> TGraph;
typedef vector<int> TIntegers;
typedef vector<double> TDoubles;
typedef queue<int> TQueue;

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	int m;
	scanf("%d%%d", &n, &m);

	int s;
	int f;
	--s;
	--f;

	TGraph g(n);
	for (int i = 0; i < m; ++i) {
		int begin;
		int end;
		int cost;
		scanf("%d%d%d", &begin, &end, &cost);
		double dCost = log(1.0 - static_cast<double>(cost)/100.0);
		g[begin].push_back( TNode(end, dCost) );
		g[end].push_back(TNode(begin, dCost));
		--begin;
		--end;
	}

	TIntegers dist(-1, n);
	TDoubles ddist(0.0, n);
	TIntegers parent(-1, n);
	dist[s] = 0;
	ddist[s] = 0.0;
	TQueue q;
	q.push(s);
	while (!q.empty()) {
		int now = q.front();
		q.pop();
		for (int i = 0; i < g[now].size(); ++i) {
			int next = g[now][i]._to;
			if (-1 == dist[next]) {
				q.push(next);
				dist[next] = dist[now] + 1;
				ddist[next] = d
			}
		}
	}

	return 0;
}