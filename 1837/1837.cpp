#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <vector>
#include <map>
#include <string>
#include <queue>

using namespace std;

typedef vector<string> TStrings;
typedef map<string, TStrings> TGraph;
typedef map<string, int> TResult;
typedef queue<string> TQueue;

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	TGraph g;

	int n;
	scanf("%d", &n);
	for (int i = 0; i < n; ++i) {
		char s1[1000];
		char s2[1000];
		char s3[1000];
		scanf("%s%s%s", s1, s2, s3);
		g[s1].push_back(s2);
		g[s1].push_back(s3);
		g[s2].push_back(s1);
		g[s2].push_back(s3);
		g[s3].push_back(s1);
		g[s3].push_back(s2);
	}

	static const string ISENBAEV = "Isenbaev";
	if (g.find(ISENBAEV) != g.end()) {
		TResult result;
		TQueue q;
		q.push(ISENBAEV);
		result[ISENBAEV] = 0;
		while (!q.empty()) {
			string top = q.front();
			q.pop();
			const TStrings& v = g[top];
			for (int i = 0; i < v.size(); ++i) {
				if (result.find(v[i]) == result.end()) {
					result[v[i]] = result[top] + 1;
					q.push(v[i]);
				}
			}
		}
	}

	for (TGraph::const_iterator toV = g.begin(); toV != g.end(); ++toV) {
		printf("%s ", toV->first.c_str());
		if (result.find(toV->first) != result.end()) {
			printf("%d\n", result[toV->first]);
		}
		else {
			printf("undefined\n");
		}
	}

	return 0;
}