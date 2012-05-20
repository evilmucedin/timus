#include <cstdio>
#include <cstring>
#include <cstdlib>

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
        if (*s >= '0' && *s <= '9') {
            return true;
        }
        ++s;
    }
    return false;
}

bool IsDelim(char ch) {
    return ch == '\r' || ch == '\n' || ch == ' ' || ch == '\t';
}

int main() {
    int n;
    scanf("%d", &n);
    char buffer[123456];
    gets(buffer);
    for (int iTest = 0; iTest < n; ++iTest) {
        gets(buffer);
        while (IsDelim(buffer[strlen(buffer) - 1]))
            buffer[strlen(buffer) - 1] = 0;
        if (!HasDigit(buffer)) {
            while (1);
        }
        if (IsDelim(*buffer)) {
            while (1);
        }
        size_t len = strlen(buffer);
        while (Rem7(buffer)) {
            size_t index1 = rand() % len;
            if (buffer[index1] == '0')
                continue;
            size_t index2 = rand() % len;
            if (buffer[index2] == '0')
                continue;
            if (index1 == index2)
                continue;
            char temp = buffer[index1];
            buffer[index1] = buffer[index2];
            buffer[index2] = temp;
        }
        printf("%s\n", buffer);
    }
    return 0;
}
