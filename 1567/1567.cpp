#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>

#include <string>

using namespace std;

int main() {
	static const string KEYBOARD[] = { "abc", "def", "ghi", "jkl", "mno", "pqr", "stu", "vwx", "yz", ".,!", " " };
	
	char line[10000];
	gets(line);
	const size_t len = strlen(line);

	int result = 0;
	for (size_t i = 0; i < len; ++i) {
		for (size_t j = 0; j < 11; ++j) {
			size_t toChar = KEYBOARD[j].find(line[i]);
			if (toChar != string::npos) {
				result += toChar  + 1;
			}
		}
	}

	printf("%d\n", result);

	return 0;
}