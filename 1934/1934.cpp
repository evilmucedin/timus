#ifndef _MSC_VER
#   pragma GCC target("sse4.2")
#   pragma GCC optimize("O3")
#endif

#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cmath>

#include <vector>
#include <queue>
#include <algorithm>
#include <random>

using namespace std;

struct TNode {
	int _to;
    double _p;

    TNode(int to, double p)
		: _to(to)
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

    bool operator<(const TCost& rhs) const {
        if (_iDist != rhs._iDist) {
            return _iDist < rhs._iDist;
        } else {
            return _fine > rhs._fine;
        }
    }
};

typedef vector<TNode> TNodes;
typedef vector<TNodes> TGraph;
typedef vector<int> TIntegers;
typedef queue<int> TQueue;
typedef vector<TCost> TCosts;

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

	int s = ReadInt() - 1;
	int f = ReadInt() - 1;

	TGraph g(n);
	for (int i = 0; i < m; ++i) {
		int begin = ReadInt() - 1;
		int end = ReadInt() - 1;
		int cost = ReadInt();

        double p = 1.0 - static_cast<double>(cost) / 100.0;
		g[begin].push_back( TNode(end, p) );
		g[end].push_back( TNode(begin, p) );
	}

    /*
    for (int i = 0; i < n; ++i) {
        shuffle(g[i].begin(), g[i].end(), default_random_engine());
    }
    */

    const TCost INF(n + 1000, -10.0);
    TCosts costs(n, INF);
    TIntegers parent(n, -1);
    parent[s] = -2;

    TQueue q;
    {
        const TCost cost(0, 1.0);
        q.push(s);
        costs[s] = cost;
    }

    while (!q.empty()) {
        const int now = q.front();
        q.pop();
        if (now == f) {
            break;
        }
        const TCost& nowCost = costs[now];
        for (TNodes::const_iterator toNode = g[now].begin(); toNode != g[now].end(); ++toNode) {
            const int next = toNode->_to;
            TCost nextCost(nowCost._iDist + 1, nowCost._fine*toNode->_p);
            if (-1 == parent[next] || ((nowCost._iDist + 1 == costs[next]._iDist) && (nowCost._fine*toNode->_p > costs[next]._fine))) {
                if (-1 == parent[next]) {
                    q.push(next);
                }
                costs[next] = nextCost;
                parent[next] = now;
            }
        }
    }

    TIntegers path;
    int now = f;
    while (-2 != now) {
        path.push_back(now);
        now = parent[now];
    }

    printf("%d %.12lf\n", path.size(), 1.0 - costs[f]._fine);
    for (TIntegers::const_reverse_iterator toPath = path.rbegin(); toPath != path.rend(); ++toPath) {
        printf("%d ", *toPath + 1);
    }
    printf("\n");

	return 0;
}