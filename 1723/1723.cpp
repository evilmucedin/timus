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

	string spell;
	cin >> spell;

	static const size_t NCOUNTS = 256;
	int counts[NCOUNTS];
	memset(counts, 0, sizeof(counts));
	for (char ch: spell)
	{
		unsigned char uch = static_cast<unsigned char>(ch);
		++counts[uch];
	}

	int max = 0;
	unsigned char chMax;
	for (size_t i = 0; i < NCOUNTS; ++i)
	{
		if (counts[i] > max)
		{
			max = counts[i];
			chMax = i;
		}
	}

	cout << chMax << endl;

	return 0;
}
