#include <iostream>

#include <cstring>

using namespace std;

int main() {
    static const int LIMIT = 37;
    int remainders[LIMIT];
    memset(remainders, 0, sizeof(remainders));
    int min = 1;
    while (cin) {
        char ch = -1;
        cin >> ch;
        if (-1 != ch) {
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
    int result = min + 1;
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
