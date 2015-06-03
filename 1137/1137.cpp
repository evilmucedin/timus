#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

using namespace std;

typedef list<int> TPath;
typedef vector<TPath> TPaths;

int main() {
	int n;
	cin >> n;
	TPaths paths(n);

	for (int i = 0; i < n; ++i) {
		int m;
		cin >> m;
		for (int j = 0; j <= m; ++j) {
			int v;
			cin >> v;
			paths[i].push_back(v);
		}
		paths[i].pop_back();
	}

	auto& result = paths[0];
	while (1 != paths.size()) {
		for (unsigned i = 1; i < paths.size();) {
			auto& path = paths[i];
			bool merged = false;
			for (auto it = result.begin(); it != result.end(); ++it) {
				auto v = find(path.begin(), path.end(), *it);
				if (v != path.end()) {
					rotate(path.begin(), v, path.end());
					result.insert(it, path.begin(), path.end());
					merged = true;
					break;
				}
			}
			if (merged) {
				paths.erase(paths.begin() + i);
			} else {
				++i;
			}
		}
	}
	result.push_back(result.front());

	cout << result.size() - 1 << " ";
	for (auto p : result) {
		cout << p << " ";
	}
	cout << endl;
	
	return 0;
}
