#include <iostream>

using namespace std;

int main()
{
    int n;
    int max = 0;
    int maxi = -1;

    cin >> n;
    int a2;
    cin >> a2;
    int a3;
    cin >> a3;
    int sum = a2 + a3;
    int a1 = 0;
    for (int i = 2; i < n; ++i)
    {
        int a;
        cin >> a;
        sum += a;
        sum -= a1;

        if (sum > max)
        {
            max = sum;
            maxi = i;
        }

        a1 = a2;
        a2 = a3;
        a3 = a;
    }

    cout << max << " " << maxi << endl;

    return 0;
}