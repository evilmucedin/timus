#include <iostream>

using namespace std;

int main()
{
    int a;
    int b;
    int c;

    cin >> a >> b >> c;

    static const int INF = 1 << 30;
    int res = INF;

    res = std::min(res, a+b+c);
    res = std::min(res, a+b-c);
    res = std::min(res, a+b*c);
    res = std::min(res, a-b+c);
    res = std::min(res, a-b-c);
    res = std::min(res, a-b*c);
    res = std::min(res, a*b+c);
    res = std::min(res, a*b-c);
    res = std::min(res, a*b*c);

    cout << res << endl;

    return 0;
}
