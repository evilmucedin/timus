#include <iostream>

#include <cstring>

using namespace std;

int main() {
    static const int LIMIT = 10 + 'Z' - 'A';
    int remainders[LIMIT];
    memset(remainders, 0, sizeof(remainders));
    int min = 2;
    while (cin) {
        char ch;
        cin >> ch;
        if (cin) {
            int value = -1;
            if (ch >= '0' && ch <= '9') {
                value = ch - '0';
            } else if (ch >= 'A' && ch <= 'Z') {
                value = ch - 'A' + 10;
            }
            if (-1 != value) {
                if (value > min) {
                    min = value;
                }
                for (size_t i = 2; i < LIMIT; ++i) {
                    remainders[i] = (remainders[i]*i + value) % (i - 1);
                }
            }
        }
    }
    int result = min;
    while (result < LIMIT && remainders[result]) {
        ++result;
    }
    if (result < LIMIT) {
        cout << result << endl;
    } else {
        cout << "No solution." << endl;
    }
    return 0;
}
