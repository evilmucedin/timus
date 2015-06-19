#include <iostream>

using namespace std;

int main()
{
	int n;
	int k;

	cin >> n >> k;
	for (int i = 0; i < n; ++i)
	{
		int b;
		int g;
		cin >> b >> g;
		k += b;
		k -= g + 2;
	}
	k -= 2;
	if (k >= 0)
	{
		cout << k << endl;
	}
	else
	{
		cout << "Big Bang!" << endl;
	}
	
	return 0;
}
