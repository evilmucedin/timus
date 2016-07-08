#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n;
    cin >> n;
    int m;
    cin >> m;
    vector<int> nums(n);
    for (auto& x: nums) {
        cin >> x;
    }
    for (int i = 0; i < 10; ++i) {
        cout << nums[(i + m) % n];
    }
    cout << endl;
    return 0;
}
