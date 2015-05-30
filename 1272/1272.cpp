#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>

using namespace std;

typedef vector<int> TIntegers;

struct TGraph {
	TIntegers _component;
	TIntegers _next;
	int _nComponents;

	TGraph(int n)
		: _component(n)
		, _next(n, -1)
		, _nComponents(n)
	{
		for (int i = 0; i < n; ++i) {
			_component[i] = i;
		}
	}

	void Unite(int i, int j) {
		int now = j;
		int prev = -1;
		while (-1 != now) {
			_component[now] = _component[i];
			prev = now;
			now = _next[now];
		}

		int oldnext = _next[i];
		_next[i] = j;
		if (-1 != oldnext) {
			_next[prev] = oldnext;
		}
		
		--_nComponents;
	}
};

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	int k;
	int m;
	scanf("%d%d%d", &n, &k, &m);

	TGraph g(n);
	for (int i = 0; i < k; ++i) {
		int begin;
		int end;
		scanf("%d%d", &begin, &end);
		--begin;
		--end;

		if (g._component[begin] != g._component[end]) {
			g.Unite(g._component[begin], g._component[end]);
		}
	}

	printf("%d\n", g._nComponents - 1);

	return 0;
}