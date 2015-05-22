#include <iostream>
#include <string>

using namespace std;

int main() {
	int n;
	cin >> n;
	int now = 0;
	int result = 0;
	for (int i = 0; i < n; ++i) {
		string s;
		cin >> s;
		int cs;
		switch (s[0]) {
		case 'A': case 'P': case 'O': case 'R':
			cs = 0;
			break;
		case 'B': case 'M': case 'S':
			cs = 1;
			break;
		default:
			cs = 2;
			break;
		}
		int diff = now - cs;
		if (diff < 0) {
			diff = -diff;
		}
		result += diff;
		now = cs;
	}
	cout << result << endl;
	return 0;
}