#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cctype>

int ReadInt() {
    int result = 0;
    char ch = getchar();
    while (!isdigit(ch))
        ch = getchar();
    while (isdigit(ch)) {
        result = 10 * result + ch - '0';
        ch = getchar();
    }
    return result;
}

int main() {
    int n = ReadInt();
    int k = ReadInt();
    int remBB = 0;
    int remDroids = 0;
    for (int i = 0; i < n; ++i) {
        int bb = ReadInt();
        if (bb > k) {
            remBB += bb - k;
        } else {
            remDroids += k - bb;
        }
    }
    printf("%d %d\n", remBB, remDroids);
    return 0;
}