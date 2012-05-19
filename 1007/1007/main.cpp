#include <cstdio>
#include <string>
#include <cassert>
#include <algorithm>

using namespace std;

bool IsCorrect(const std::string& s, int n) {
    assert(s.length() == n);
    int res = 0;
    for (size_t i = 0; i < s.length(); ++i)
        if (s[i] == '1')
            res += i + 1;
    return 0 == (res % (n + 1));
}

int main() {
    freopen("input.txt", "r", stdin);
    
    int n;
    scanf("%d", &n);
    char buffer[10000];
    while (1 == scanf("%s", buffer)) {
        string s = buffer;

        int mod;
        {
            int res = 0;
            for (size_t i = 0; i < s.length(); ++i)
                if (s[i] == '1')
                    res += i + 1;
            mod = res % (n + 1);
        }
        
        string res;
        if (s.length() == n) {
            if (IsCorrect(s, n)) {
                res = s;
            } else {
                for (size_t i = max(mod - 2, 0); i < min(s.length(), (size_t)mod + 2); ++i) {
                    string cand = s;
                    cand[i] = '0';
                    if (IsCorrect(cand, n))
                        res = cand;
                }
            }
        } else if (s.length() == n + 1) {
            for (size_t i = max(mod - 10, 1); i < min(s.length(), (size_t)mod + 10); ++i) {
                string cand = s.substr(0, i - 1) + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    // printf("%d %d\n", i, mod);
                    res = cand;
                }
            }            
            for (size_t i = max(n - mod - 10, 1); i < min(s.length(), n - (size_t)mod + 10); ++i) {
                string cand = s.substr(0, i - 1) + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    // printf("%d %d\n", i, mod);
                    res = cand;
                }
            }            
        } else if (s.length() == n - 1) {
            for (size_t i = 0; i <= 0; ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = s.length(); i <= s.length(); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(mod - 5, 0); i <= min(s.length(), (size_t)mod + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(n - mod - 5, 0); i <= min(s.length(), n - (size_t)mod + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(mod/2 - 5, 0); i <= min(s.length(), (size_t)mod/2 + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(n - mod/2 - 5, 0); i <= min(s.length(), n - (size_t)mod/2 + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(mod*2 - 5, 0); i <= min(s.length(), (size_t)mod*2 + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
            for (size_t i = max(n - mod*2 - 5, 0); i <= min(s.length(), n - (size_t)mod*2 + 5); ++i) {
                string cand = s.substr(0, i) + "0" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
                cand = s.substr(0, i) + "1" + s.substr(i, s.length());
                if (IsCorrect(cand, n)) {
                    res = cand;
                }
            }                        
        }
        printf("%s\n", res.c_str());
    }
    
    return 0;
}