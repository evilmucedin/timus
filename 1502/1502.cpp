#include <iostream>

using namespace std;

int main() {
    int n;
    cin >> n;
    long long int result = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = i; j <= n; ++j) {
            result += i + j;
        }
    }
    cout << result << endl;
    return 0;
}
