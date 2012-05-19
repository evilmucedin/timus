#include <cstdio>
#include <string>
#include <cassert>
#include <algorithm>
#include <vector>

using namespace std;

typedef vector<int> TIntVector;

bool IsGood(int num, int n) {
    while (num < 0)
        num += n + 1;
    return (num % (n + 1)) == 0;
}

int main() {
    freopen("input.txt", "r", stdin);
    
    int n;
    scanf("%d", &n);
    char buffer[10000];
    while (1 == scanf("%s", buffer)) {
        string s = buffer;

        int mod;
        TIntVector counts(s.length() + 1);
        {
            int res = 0;
            for (size_t i = 0; i < s.length(); ++i)
                if (s[i] == '1')
                    res += i + 1;
            mod = res % (n + 1);
            counts[s.length()] = 0;
            for (int i = (int)s.length() - 1; i >= 0; --i) {
                counts[i] = counts[i + 1];
                if (s[i] == '1')
                    ++counts[i];
            }
        }
        
        string res;
        if (s.length() == n) {
            if (IsGood(mod, n)) {
                res = s;
            } else {
                for (int i = 0; i < n; ++i) {
                    if (s[i] == '1') {
                        if ( IsGood(mod - i - 1, n) ) {
                            res = s;
                            res[i] = '0';
                        }
                    }
                }
            }
        } else if (s.length() == n + 1) {
            for (int i = 0; i < s.length(); ++i) {
                if (s[i] == '1') {
                    if (IsGood(mod - i - 1 - counts[i + 1], n))
                        res = s.substr(0, i) + s.substr(i + 1, s.length());
                } else {
                    if (IsGood(mod - counts[i + 1], n))
                        res = s.substr(0, i) + s.substr(i + 1, s.length());                    
                }
            }
        } else if (s.length() == n - 1) {
            for (int i = 0; i <= s.length(); ++i) {
                if (IsGood(mod + i + 1 + counts[i], n))
                    res = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsGood(mod + counts[i], n))
                    res = s.substr(0, i) + "0" + s.substr(i, s.length());
            }
        }
        printf("%s\n", res.c_str());
    }
    
    return 0;
}