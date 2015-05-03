#include <cstdio>
#include <cstring>
#include <cstdlib>

bool IsDigit(char ch) {
    return ch >= '0' && ch <= '9';
}

int Rem7(const char* buffer) {
    int rem = 0;
    while (*buffer) {
        rem = (10*rem + (((int)*buffer) - '0')) % 7;
        ++buffer;
    }
    return rem;
}

bool HasDigit(const char* s) {
    while (*s) {
        if (IsDigit(*s))
            return true;
        ++s;
    }
    return false;
}

size_t IndexOf(const char* buffer, char ch) {
    size_t index = 0;
    while (buffer[index] != ch)
        ++index;
    return index;
}

int main() {
    int n;
    scanf("%d", &n);
    char buffer[1234];
    gets(buffer);
    for (int iTest = 0; iTest < n; ++iTest) {
        gets(buffer);

        size_t index[4];
        for (size_t i = 0; i < 4; ++i)
            index[i] = IndexOf(buffer, i + 1 + '0');

        static const int PERM[24][4] = { {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2}, {0, 3, 2, 1},
                                      {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0},
                                      {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0},
                                      {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0} };
        char copy[1234];
        size_t i = 0;
        for (i = 0; i < 24; ++i) {
            strcpy(copy, buffer);
            for (size_t j = 0; j < 4; ++j)
                copy[ index[ PERM[i][j] ] ] = buffer[index[j]];
            if (!Rem7(copy))
                break;
        }

        if (i == 24) {
            size_t len = strlen(copy);
            while (copy[0] == '0' || Rem7(copy)) {
                size_t index1 = rand() % len;
                size_t index2 = rand() % len;
                char temp = copy[index1];
                copy[index1] = copy[index2];
                copy[index2] = temp;
            }
        }

        printf("%s\n", copy);
    }
    return 0;
}
