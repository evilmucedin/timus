#include <memory.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
#ifndef ONLINE_JUDGE
	ifstream fIn("input.txt");
	cin.rdbuf(fIn.rdbuf());
#endif

	string s;
	cin >> s;

	int counts[256];
	memset(counts, 0, sizeof(counts));
	for (char ch: s)
	{
		unsigned char uch = static_cast<unsigned char>(ch);
		++counts[uch];
	}

	int max = 0;
	unsigned char iMax;
	for (int i = 0; i < 256; ++i)
	{
		if (counts[i] > max)
		{
			max = counts[i];
			iMax = i;
		}
	}

	cout << iMax << endl;

	return 0;
}
