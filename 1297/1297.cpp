#include <string>
#include <iostream>
#include <vector>

using namespace std;

typedef vector<bool> TBoolVector;
typedef vector<TBoolVector> TBoolMatrix;

int main() {
    string s;
    cin >> s;

    TBoolVector dummy(s.length() + 1);
    TBoolMatrix matrix(s.length() + 1, dummy);

    for (size_t i = 0; i < s.length(); ++i) {
        matrix[i][i] = true;
    }
    for (size_t i = 0; i < s.length(); ++i) {
        matrix[i][i + 1] = true;
    }
    for (int len = 2; len <= s.length(); ++len) {
        for (int i = 0; i <= s.length() && i + len <= s.length(); ++i) {
            if (matrix[i + 1][i + len - 1] && s[i] == s[i + len - 1]) {
                matrix[i][i + len] = true;
            }
        }
    }

    int max = 0;
    int maxi = 0;
    int maxj = 0;
    for (int i = 0; i <= s.length(); ++i) {
        for (int j = i; j <= s.length(); ++j) {
            if (matrix[i][j]) {
                if (j - i > max) {
                    max = j - i;
                    maxi = i;
                    maxj = j;
                }
            }
        }
    }

    cout << s.substr(maxi, max) << endl;

    return 0;
}