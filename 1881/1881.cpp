#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstring>

int main() {
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
#endif

    int h;
    int w;
    int n;
    scanf("%d%d%d", &h, &w, &n);

    int nPages = 0;
    int nLine = h;
    int nChars = 0;

    for (int i = 0; i < n; ++i) {
        char buffer[1000];
        scanf("%s", buffer);
        int len = strlen(buffer);
        if (nChars && nChars + len + 1 <= w) {
            nChars += len + 1;
        } else {
            if (nLine == h) {
                ++nPages;
                nLine = 0;
            }
            ++nLine;
            nChars = len;
        }
    }

    printf("%d\n", nPages);
    
    return 0;
}