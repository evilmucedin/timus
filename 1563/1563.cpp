#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <set>
#include <string>

using namespace std;

typedef set<string> TStringSet;

int main()
{
	freopen("input.txt", "r", stdin);

	int n;
	scanf("%d", &n);

	TStringSet visited;

	int skipped = 0;

	char buffer[1000];
	gets(buffer);
	while (gets(buffer))
	{
		string sBuffer(buffer);
		if (visited.end() != visited.find(sBuffer))
		{
			++skipped;
		}
		else
		{
			visited.insert(sBuffer);
		}
	}

	printf("%d\n", skipped);

	return 0;
}