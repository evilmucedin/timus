#include <iostream>

#include <string>
#include <algorithm>

using namespace std;

int main()
{
    int t;
    cin >> t;

    for (int iTest = 0; iTest < t; ++iTest)
    {
        long long int n;
        long long int k;
        cin >> n >> k;

        long long int groupSize = n/k;
        long long int rem = n%k;

        long long int answer = ((n - (groupSize + 1)*rem)*(n - groupSize) + (rem*(groupSize + 1))*(n - groupSize - 1))/2; 

        cout << answer << endl;
    }

    return 0;
}