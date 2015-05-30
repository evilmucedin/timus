#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>
#include <queue>

using namespace std;

typedef vector<int> TIntegers;

typedef vector<char> TBoolLine;
typedef vector<TBoolLine> TBoolMatrix;

typedef queue<int> TIntQueue;

char ReadNext() {
	char ch = getchar();
	while (ch != '0' && ch != '1')
		ch = getchar();
	return ch;
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif
	
	int n;
	scanf("%d", &n);
	int a, d;
	scanf("%d%d", &d, &a);

	TBoolMatrix graph(n, TBoolLine(n));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			graph[i][j] = ReadNext();
		}
	}

	long long int cost = 0;
	TIntegers components(n, -1);
	int nComponents = 0;
	for (int i = 0; i < n; ++i) {
		if (-1 == components[i]) {
			TIntQueue q;
			q.push(i);
			components[i] = nComponents;
			while (!q.empty()) {
				int v = q.front();
				q.pop();
				for (int j = 0; j < n; ++j) {
					if ('1' == graph[v][j]) {
						if (-1 == components[j]) {
							q.push(j);
							components[j] = nComponents;
							graph[v][j] = '0';
							graph[j][v] = '0';
						} else {
							cost += d;
							graph[v][j] = 'd';
							graph[j][v] = 'd';
						}
					}
				}
			}
			++nComponents;
		}
	}

	int iComponent = 1;
	for (int i = 1; i < n; ++i) {
		if (components[i] == iComponent) {
			graph[0][i] = 'a';
			graph[i][0] = 'a';
			++iComponent;
			cost += a;
		}
	}

	printf("%lld\n", cost);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%c", graph[i][j]);
		}
		printf("\n");
	}

	return 0;
}