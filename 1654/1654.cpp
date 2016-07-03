#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <stack>
#include <vector>

using namespace std;

typedef vector<char> TCharVector;

int main()
{
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
#endif
	char ch;

	TCharVector v;
	while (-1 != (ch = getchar())) 
	{
		if (ch >= 'a' && ch <= 'z')
		{
			if (v.empty() || v.back() != ch)
			{
				v.push_back(ch);
			}
			else
			{
				v.pop_back();
			}
		}
	}

	for (int i = 0; i < v.size(); ++i)
	{
		printf("%c", v[i]);
	}
	printf("\n");

	return 0;
}