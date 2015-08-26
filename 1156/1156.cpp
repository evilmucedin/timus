#include <iostream>
#include <fstream>

#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

typedef vector<bool> TBoolVector;
typedef vector<TBoolVector> TMatrix;

typedef vector<int> TIntVector;

int main() 
{
#ifndef ONLINE_JUDGE
	ifstream fIn("input.txt");
	cin.rdbuf( fIn.rdbuf() );
#endif

	int n;
	int m;
	cin >> n >> m;

	TBoolVector dummy(2*n, false);
	TMatrix g(2*n, dummy);
	for (int i = 0; i < m; ++i) 
	{
		int a;
		int b;
		cin >> a >> b;
		--a;
		--b;
		g[a][b] = true;
		g[b][a] = true;
	}

	TIntVector s;
	s.push_back(0);
	do 
	{
		int i = s.back() + 1;
		for (; i < 2 * n; ++i)
		{
			int j = 0;
			while (j < s.size() && !g[s[j]][i])
			{
				++j;
			}
			if (j == s.size())
			{
				break;
			}
		}
		if (i == 2 * n)
		{
			s.pop_back();
		}
		else
		{
			s.push_back(i);
		}
	} while (!s.empty() && s.size() != n);
	
	if (s.size() == n)
	{
		for (int i = 0; i < 2 * n; ++i)
		{
			if (find(s.begin(), s.end(), i) != s.end())
			{
				cout << i + 1 << " ";
			}
		}
		cout << endl;
		for (int i = 0; i < 2 * n; ++i)
		{
			if (find(s.begin(), s.end(), i) == s.end())
			{
				cout << i + 1 << " ";
			}
		}
		cout << endl;
	}
	else
	{
		cout << "IMPOSSIBLE\n";
	}
	
	return 0;
}