#include <iostream>

#include <map>

using namespace std;

typedef map<int, int> TIntIntMap;

int main()
{
    TIntIntMap counts;
    for (int i = 0; i < 3; ++i)
    {
        int n;
        cin >> n;
        for (int j = 0; j < n; ++j)
        {
            int ev;
            cin >> ev;
            ++counts[ev];
        }
    }

    int count = 0;
    for (TIntIntMap::const_iterator toCount = counts.begin(); toCount != counts.end(); ++toCount)
    {
        if (toCount->second == 3)
            ++count;
    }

    cout << count << endl;

    return 0;
}