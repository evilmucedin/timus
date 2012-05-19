#include <cstdio>
#include <cmath>

#pragma comment(linker, "/stack:10000000")

void run() {
    double v;
    if (1 == scanf("%lf", &v)) {
        run();
        printf("%.20lf\n", sqrt(v));
    }
}

int main()
{
    run();
    return 0;
}

