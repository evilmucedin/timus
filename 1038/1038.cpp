#include <cstdio>

bool IsSmall(char ch) {
    return (ch >= 'a') && (ch <= 'z');
}

bool IsBig(char ch) {
    return (ch >= 'A') && (ch <= 'Z');
}

bool IsLetter(char ch) {
    return IsSmall(ch) || IsBig(ch);
}

bool IsEnd(char ch) {
    return ch == '.' || ch == '?' || ch == '!';
}

int main() {
    int mistakes = 0;
    char ch;
    bool inWord = false;
    bool sentenceBegin = true;
    while (1 == scanf("%c", &ch)) {
        if (IsLetter(ch)) {
            if (inWord)
                if (IsBig(ch))
                    ++mistakes;
            if (sentenceBegin)
                if (IsSmall(ch))
                    ++mistakes;
            inWord = true;
            sentenceBegin = false;
        } else {
            inWord = false;
            if (IsEnd(ch))
                sentenceBegin = true;
        }
    }
    printf("%d\n", mistakes);
    return 0;
}
