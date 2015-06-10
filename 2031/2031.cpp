#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>

#include <string>

using namespace std;

int main()
{
    int n;
    scanf("%d", &n);
    if (n > 4)
    {
        printf("Glupenky Pierre\n");
    }
    else
    {
        static const string data[] = {"16", "06", "68", "88"};
        for (int i = 0; i < n; ++i)
        {
            printf("%s ", data[i].c_str());
        }
        printf("\n");
    }
    return 0;
}