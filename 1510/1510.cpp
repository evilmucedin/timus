#include <cstdio>

int main() {
    int n;
    scanf("%d", &n);

    int candidate = 0;
    int counter = 0;
    for (int i = 0; i < n; ++i) {
        int x;
        scanf("%d", &x);
        if (x == candidate) {
            ++counter;
        } else {
            --counter;
        }
        if (counter < 0) {
            candidate = x;
            counter = 1;
        }
    }

    printf("%d\n", candidate);

    return 0;
}
