#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>
#include <memory>

using namespace std;

class Exception : public std::exception {
private:
    const char* _s;

public:
    Exception(const char* s)
        : _s(s)
    {
    }

    const char* what() const throw() {
        return _s;
    }
};

struct IParser {
    virtual bool Parse(char*& s) = 0;
};

typedef vector<bool> TBoolVector;
typedef vector<TBoolVector> TBoolMatrix;

struct ITemplateParser : public IParser {
    virtual void FillMask(TBoolVector* vct) {

    }

    virtual bool IsStar() {
        return false;
    }
};

typedef shared_ptr<ITemplateParser> PTemplateParser;

struct TCharacterParser : public ITemplateParser {
    char _ch;

    virtual void FillMask(TBoolVector* vct) {
        (*vct)[UChar()] = true;
    }

    bool Parse(char*& s) {
        if (*s == '\'') {
            if (*(s + 1) == '\'') {
                s += 2;
                _ch = '\'';
                return true;
            }
        } else {
            unsigned char uc = *s;
            if (uc >= 32) {
                _ch = *s;
                ++s;
                return true;
            }
        }
        return false;
    }

    unsigned char UChar() {
        return static_cast<unsigned char>(_ch);
    }
};

struct TStringParser : public IParser {
    vector<TCharacterParser> _chars;

    bool Parse(char*& s) {
        char* backup = s;
        while (*s) {
            TCharacterParser cp;
            if (cp.Parse(s)) {
                _chars.push_back(cp);
            } else {
                break;
            }
        }
        if (!*s) {
            s = backup;
            return false;
        }
        return true;
    }
};

struct TPercentParser : public ITemplateParser {
    virtual void FillMask(TBoolVector* vct) {
    }

    virtual bool IsStar() {
        return true;
    }

    bool Parse(char*& s) {
        if (*s == '%') {
            ++s;
            return true;
        }
        return false;
    }
};

struct TUnderscoreParser : public ITemplateParser {
    virtual void FillMask(TBoolVector* vct) {
        for (size_t i = 0; i < vct->size(); ++i) {
            (*vct)[i] = true;
        }
    }

    bool Parse(char*& s) {
        if (*s == '_') {
            ++s;
            return true;
        }
        return false;
    }
};

struct TRangeParser : public ITemplateParser {
    TCharacterParser _c1;
    TCharacterParser _c2;

    virtual void FillMask(TBoolVector* vct) {
        for (size_t i = _c1.UChar(); i <= _c2.UChar(); ++i) {
            (*vct)[i] = true;
        }
    }

    bool Parse(char*& s) {
        char* backup = s;
        if (_c1.Parse(s)) {
            if (*s == '-') {
                ++s;
                if (_c2.Parse(s)) {
                    return true;
                }
            }
        }
        s = backup;
        return false;
    }
};

struct TSetParser : ITemplateParser {
    bool _not;
    vector<PTemplateParser> _parsers;

    virtual void FillMask(TBoolVector* vct) {
        for (size_t i = 0; i < _parsers.size(); ++i) {
            _parsers[i]->FillMask(vct);
        }
        if (_not) {
            for (size_t i = 0; i < vct->size(); ++i) {
                (*vct)[i] = !(*vct)[i];
            }
        }
    }

    bool Parse(char*& s) {
        _not = false;
        char* backup = s;
        if (*s == '[') {
            ++s;

            if (*s == '^') {
                _not = true;
                ++s;
            }

            while (*s && *s != ']') {
                PTemplateParser parser = make_shared<TRangeParser>();
                if (parser->Parse(s)) {
                    _parsers.push_back(parser);
                    continue;
                }

                parser = make_shared<TCharacterParser>();
                if (parser->Parse(s)) {
                    _parsers.push_back(parser);
                    continue;
                }

                s = backup;
                return false;
            }
            if (*s == ']') {
                ++s;
                return true;
            }
        }
        s = backup;
        return false;
    }
};

struct TTemplateItem {
    char _isStar;
    TBoolVector _mask;

    TTemplateItem()
        : _mask(256)
    {
    }
};

typedef vector<TTemplateItem> TTemplateItems;

struct TTemplateParser {
    typedef shared_ptr<ITemplateParser> PTemplateParser;
    vector<PTemplateParser> _parsers;

    void ToItems(TTemplateItems* items) {
        items->clear();
        items->resize(_parsers.size());
        for (size_t i = 0; i < _parsers.size(); ++i) {
            (*items)[i]._isStar = _parsers[i]->IsStar();
            if (!(*items)[i]._isStar) {
                _parsers[i]->FillMask(&((*items)[i]._mask));
            }
        }
    }

    bool Parse(char*& s) {
        char* backup = s;
        while (*s) {
            PTemplateParser pp = make_shared<TPercentParser>();
            if (pp->Parse(s)) {
                _parsers.push_back(pp);
                continue;
            }
            pp = make_shared<TUnderscoreParser>();
            if (pp->Parse(s)) {
                _parsers.push_back(pp);
                continue;
            }
            pp = make_shared<TSetParser>();
            if (pp->Parse(s)) {
                _parsers.push_back(pp);
                continue;
            }
            pp = make_shared<TCharacterParser>();
            if (pp->Parse(s)) {
                _parsers.push_back(pp);
                continue;
            }
            break;
        }
        if (!*s) {
            s = backup;
            return false;
        }
        return true;
    }
};

struct TLineParser : public IParser {
    TStringParser _string;
    TTemplateParser _template;

    bool Parse(char*& s) {
        if (*s != '\'') {
            return false;
        }
        ++s;
        if (!_string.Parse(s)) {
            return false;
        }
        if (0 != memcmp(s, "\' like \'", 8)) {
            return false;
        }
        s += 8;
        if (!_template.Parse(s)) {
            return false;
        }
        if (*s != '\'') {
            return false;
        }
        return true;
    }

    static bool Get(const TBoolMatrix& m, int i, int j) {
        if (i < 0) {
            return false;
        }
        if (j < 0) {
            return false;
        }
        return m[i][j];
    }

    bool Like() {
        TTemplateItems templateItems;
        _template.ToItems(&templateItems);
        TBoolVector dummy(templateItems.size() + 1);
        TBoolMatrix matrix(_string._chars.size() + 1, dummy);
        matrix[0][0] = true;
        for (int i = 0; i <= _string._chars.size(); ++i) {
            for (int j = 1; j <= templateItems.size(); ++j) {
                if (templateItems[j - 1]._isStar) {
                    matrix[i][j] = Get(matrix, i - 1, j) || Get(matrix, i - 1, j - 1) || Get(matrix, i, j - 1);
                } else {
                    if (i >= 1) {
                        matrix[i][j] = templateItems[j - 1]._mask[_string._chars[i - 1].UChar()] && Get(matrix, i - 1, j - 1);
                    } else {
                        matrix[i][j] = Get(matrix, i - 1, j - 1);
                    }
                }
            }
        }
        return matrix[_string._chars.size()][templateItems.size()];
    }
};

int main() {
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
#endif

    int n;
    scanf("%d", &n);
    char line[10000];
    for (int iTest = 0; iTest < n; ++iTest) {
        do {
            gets(line);
        } while (0 == strlen(line));
        TLineParser parser;
        char* pLine = line;
        if (!parser.Parse(pLine)) {
            throw Exception("can not parse");
        }
        printf( parser.Like() ? "YES\n" : "NO\n" );
    }

    return 0;
}