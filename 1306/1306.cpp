#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <queue>
#include <functional>

using namespace std;

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif

	int n;
	scanf("%d", &n);
	priority_queue<int, vector<int>, greater<int>> q;
	
	int med = n / 2 + 1;
	for (int i = 0; i < med; ++i)
	{
		int temp;
		scanf("%d", &temp);
		q.push(temp);
	}

	for (int i = med; i < n; ++i)
	{
		int temp;
		scanf("%d", &temp);
		q.push(temp);
		q.pop();
	}

	double result;
	if (n & 1)
	{
		result = q.top();
	}
	else
	{
		result = q.top();
		q.pop();
		result += q.top();
		result /= 2;
	}

	printf("%.1lf\n", result);
	
	return 0;
}