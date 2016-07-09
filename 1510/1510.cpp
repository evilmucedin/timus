#include <cstdio>

#include <vector>
#include <algorithm>

using namespace std;

int main() {
    int n;
    scanf("%d", &n);

    vector<int> v(n);
    for (auto& x: v) {
        scanf("%d", &x);
    }

    sort(v.begin(), v.end());

    printf("%d\n", v[n/2]);

    return 0;
}
