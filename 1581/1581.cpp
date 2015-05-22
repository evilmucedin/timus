#include <iostream>
#include <vector>

using namespace std;

typedef vector<int> Integers;

int main() {
	int n;
	cin >> n;
	Integers numbers(n);
	for (int i = 0; i < n; ++i) {
		cin >> numbers[i];
	}
	int count = 0;
	int num = -1;
	for (int i = 0; i < n; ++i) {
		if (numbers[i] == num) {
			++count;
		} else {
			if (count) {
				cout << count << " " << num << " ";
			}
			count = 1;
			num = numbers[i];
		}
	}
	if (count) {
		cout << count << " " << num;
	}
	return 0;
}