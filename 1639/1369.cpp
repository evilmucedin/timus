#include <iostream>

using namespace std;

int main() {
	int m;
	int n;
	cin >> m >> n;
	int k = m*n - 1;
	cout << ((k & 1) ? "[:=[first]" : "[second]=:]") << endl;
	return 0;
}