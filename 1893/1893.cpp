#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>

using namespace std;

enum TType {tWindow, tAsile, tUnknown};

int main()
{
	string s;
	cin >> s;
	char letter = s[s.length() - 1];
	int row;
	sscanf(s.c_str(), "%d", &row);

	TType result = tUnknown;
	if (row <= 2) {
		switch (letter) {
		case 'A':
			result = tWindow;
			break;
		case 'B':
			result = tAsile;
			break;
		case 'C':
			result = tAsile;
			break;
		case 'D':
			result = tWindow;
			break;
		}
	}
	else if (row <= 20) {
		if (letter == 'A' || letter == 'F')
			result = tWindow;
		else
			result = tAsile;
	}
	else {
		if (letter == 'A' || letter == 'K')
			result = tWindow;
		else if (letter == 'C' || letter == 'D' || letter == 'G' || letter == 'H')
			result = tAsile;
	}

	switch (result) {
	case tWindow:
		cout << "window";
		break;
	case tAsile:
		cout << "aisle";
		break;
	case tUnknown:
		cout << "neither";
		break;
	}
	cout << endl;

	return 0;
}