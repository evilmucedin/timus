#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

int main() {
    int n;
    int k;
    scanf("%d%d", &n, &k);
    int remBB = 0;
    int remDroids = 0;
    for (int i = 0; i < n; ++i) {
        int bb;
        scanf("%d", &bb);
        if (bb > k) {
            remBB += bb - k;
        } else {
            remDroids += k - bb;
        }
    }
    printf("%d %d\n", remBB, remDroids);
    return 0;
}