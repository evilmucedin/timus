#include <iostream>

using namespace std;

int main()
{
    int a;
    int b;
    cin >> a >> b;

    cout << (((a & 1) == 0 || (b & 1) == 1) ? "yes" : "no") << endl;

    return 0;
}