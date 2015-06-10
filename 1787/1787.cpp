#include <iostream>

using namespace std;

int main()
{
    int k;
    cin >> k;
    int n;
    cin >> n;
    int res = 0;
    for (int i = 0; i < n; ++i)
    {
        int delta;
        cin >> delta;
        res += delta;
        if (res < k)
            res = 0;
        else
            res -= k;
    }
    cout << res << endl;
    return 0;
}