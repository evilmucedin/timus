#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cctype>

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

	bool operator<(const TEdge& rhs) const {
		return _cost < rhs._cost;
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
	int count = 1;
	while (!stack.empty()) {
		int now = stack.back();
		stack.pop_back();
		for (TEdges::const_iterator toEdge = g[now].begin(); toEdge != g[now].end(); ++toEdge) {
			if (toEdge->_cost > ml) {
				break;
			}
			if (!visited[toEdge->_to]) {
				visited[toEdge->_to] = true;
				stack.push_back(toEdge->_to);
				++count;
			}
		}
	}
	return g.size() == count;
}

int ReadInt() {
	int result = 0;
	char ch = getchar();
	while (!isdigit(ch))
		ch = getchar();
	while (isdigit(ch)) {
		result = 10 * result + ch - '0';
		ch = getchar();
	}
	return result;
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n = ReadInt();
	int m = ReadInt();

	int mn = 100000000;
	int mx = 0;
	TGraph g(n);
	vector< tuple<int, int, int> > input(m);
	for (int i = 0; i < m; ++i) {
		int a = ReadInt();
		int b = ReadInt();
		int len = ReadInt();
		input[i] = make_tuple(a, b, len);
		--a;
		--b;
		g[a].push_back(TEdge(b, len));
		g[b].push_back(TEdge(a, len));

		mn = min(mn, len);
		mx = max(mx, len);
	}

	for (int i = 0; i < n; ++i) {
		sort(g[i].begin(), g[i].end());
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
		while ((l < r) && !Good(g, l)) {
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
	for (const auto& p : result) {
		printf("%d %d\n", p.first, p.second);
	}

	return 0;
}