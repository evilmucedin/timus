#include <iostream>
#include <fstream>

#include <vector>

using namespace std;

typedef vector<int> TIntVector;
typedef vector<TIntVector> TIntMatrix;

int main() {
#ifndef ONLINE_JUDGE
	ifstream fIn("input.txt");
	cin.set_rdbuf(fIn.rdbuf());
#endif

	TIntMatrix m(5, TIntVector(5));
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			cin >> m[i][j];
		}
	}

	static const int perms[][5] = { {0, 2, 1, 3, 4}, {0, 2, 3, 1, 4}, {0, 1, 2, 3, 4}, {0, 3, 2, 1, 4} };

	int min = 123456;
	int minIndex = -1;
	for (int i = 0; i < 4; ++i) {
		int dist = 0;
		const int* p = perms[i];
		for (int j = 0; j < 4; ++j) {
			dist += m[p[j]][p[j + 1]];
		}
		if (dist < min) {
			min = dist;
			minIndex = i;
		}
	}

	cout << min << endl;
	for (int i = 0; i < 5; ++i) {
		cout << perms[minIndex][i] + 1 << " ";
	}
	cout << endl;

	return 0;
}