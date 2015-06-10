#include <iostream>
#include <algorithm>

using namespace std;

int main()
{
    int x;
    int y;
    int c;
    cin >> x >> y >> c;

    if (x + y < c)
    {
        cout << "Impossible";
    }
    else
    {
        x = min(x, c);
        cout << x << " " << c - x << endl;
    }
    
    return 0;
}