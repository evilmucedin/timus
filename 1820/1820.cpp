#include <iostream>

using namespace std;

int main()
{
    int n;
    int k;
    cin >> n >> k;

    n *= 2;
    int ans;
    if (n % k)
        ans = n/k + 1;
    else
        ans = n/k;
    if (ans < 2)
        ans = 2;
    cout << ans << endl;

    return 0;
}