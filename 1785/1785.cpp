#include <iostream>

using namespace std;

int main()
{
    int n;
    cin >> n;
    if (n <= 4)
        cout << "few";
    else if (n <= 9)
        cout << "several";
    else if (n <= 19)
        cout << "pack";
    else if (n <= 49)
        cout << "lots";
    else if (n <= 99)
        cout << "horde";
    else if (n <= 249)
        cout << "throng";
    else if (n <= 499)
        cout << "swarm";
    else if (n <= 999)
        cout << "zounds";
    else
        cout << "legion";
    cout << endl;
    return 0;
}