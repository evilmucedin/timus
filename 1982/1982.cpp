#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

typedef vector<int> TIntegers;
typedef vector<TIntegers> TGraph;

struct TEdge {
	int _begin;
	int _end;
	int _cost;

	bool operator<(const TEdge& e) const {
		return _cost < e._cost;
	}
};

typedef vector<TEdge> TEdges;

using namespace std;

int main() {
#ifndef ONLINE_JUDGE
	ifstream fIn("input.txt");
	cin.rdbuf(fIn.rdbuf());
#endif

	int n;
	cin >> n;
	int k;
	cin >> k;

	TIntegers components(n);
	for (int i = 0; i < n; ++i) {
		components[i] = -1;
	}

	for (int i = 0; i < k; ++i) {
		int v;
		cin >> v;
		--v;
		components[v] = i;
	}

	TEdges edges;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			TEdge e;
			cin >> e._cost;
			e._begin = i;
			e._end = j;
			edges.push_back(e);
		}
	}
	sort(edges.begin(), edges.end());

	int result = 0;
	for (int i = 0; i < edges.size(); ++i) {
		const TEdge& e = edges[i];
		if (components[e._begin] == -1 && components[e._end] != -1) {
			result += e._cost;
			components[e._begin] = components[e._end];
		}
	}

	cout << result << endl;

	return 0;
}