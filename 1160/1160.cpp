#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

struct TEdge {
	int _to;
	int _cost;

	TEdge() {

	}

	TEdge(int to, int cost)
		: _to(to)
		, _cost(cost)
	{

	}
};

typedef vector<TEdge> TEdges;
typedef vector<TEdges> TGraph;

bool Good(const TGraph& g, int ml) {
	vector<bool> visited(g.size());
	vector<int> stack;
	stack.reserve(g.size());
	visited[0] = true;
	stack.push_back(0);
	while (!stack.empty()) {
		int now = stack.back();
		stack.pop_back();
		for (int i = 0; i < g[now].size(); ++i) {
			if (!visited[g[now][i]._to] && g[now][i]._cost <= ml) {
				visited[g[now][i]._to] = true;
				stack.push_back(g[now][i]._to);
			}
		}
	}
	for (int i = 0; i < g.size(); ++i) {
		if (!visited[i]) {
			return false;
		}
	}
	return true;
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	int m;
	scanf("%d%d", &n, &m);

	int mn = 100000000;
	int mx = 0;
	TGraph g(n);
	vector< tuple<int, int, int> > input;
	for (int i = 0; i < m; ++i) {
		int a;
		int b;
		int len;
		scanf("%d%d%d", &a, &b, &len);
		input.push_back(make_tuple(a, b, len));
		--a;
		--b;
		g[a].push_back(TEdge(b, len));
		g[b].push_back(TEdge(a, len));

		mn = min(mn, len);
		mx = max(mx, len);
	}

	int l = mn - 1;
	int r = mx;

	while (l + 1 < r) {
		int mid = (l + r) >> 1;
		if (Good(g, mid)) {
			r = mid;
		}
		else {
			l = mid;
		}
	}

	++l;
	if (l != r) {
		while (!Good(g, l)) {
			++l;
		}
	}

	printf("%d\n", l);

	vector< pair<int, int> > result;
	result.reserve(input.size());
	for (int i = 0; i < input.size(); ++i) {
		if (get<2>(input[i]) <= l) {
			result.push_back(make_pair(get<0>(input[i]), get<1>(input[i])));
		}
	}

	printf("%d\n", result.size());
	for (auto p : result) {
		printf("%d %d\n", p.first, p.second);
	}

	return 0;
}