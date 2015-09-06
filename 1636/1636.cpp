#include <iostream>

using namespace std;

int main()
{
	int t1;
	int t2;
	cin >> t1 >> t2;
	int sum = 0;
	for (int i = 0; i < 10; ++i)
	{
		int attempts;
		cin >> attempts;
		sum += attempts;
	}
	if (t2 >= t1 + 20*sum)
	{
		cout << "No chance." << endl;
	}
	else
	{
		cout << "Dirty debug :(" << endl;
	}
	return 0;
}
