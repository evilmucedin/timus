#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#ifdef __linux__
#   include <x86intrin.h>

typedef int HANDLE;
typedef int DWORD;
typedef void* LPVOID;
#   define WINAPI

#else
#   include <immintrin.h>
#   include <intrin.h>
#   include <windows.h>
#endif

#ifndef _MSC_VER
#   pragma GCC target("sse4.2")
#   pragma GCC optimize("O3")
#endif

#include <cstdio>
#include <cstring>
#include <cmath>

#include <vector>
#include <algorithm>

using namespace std;

typedef vector<int> TIntVector;

template<typename T>
T Sqr(T x)
{
    return x*x;
}

void Output(const TIntVector& vct)
{
    char buffer[100];
    for (size_t i = 0; i < vct.size(); ++i)
    {
        int num = vct[i] + 1;
        char* pBuffer = buffer;
        while (num)
        {
            *pBuffer = (num % 10) + '0';
            num /= 10;
            ++pBuffer;
        }
        *pBuffer = 0;
        std::reverse(buffer, pBuffer);
        *pBuffer = ' ';
        ++pBuffer;
        *pBuffer = 0;
        fputs(buffer, stdout);
    }
    fputs("\n", stdout);
}

#ifndef _MSC_VER
#   define GET_ITEM(var, index) (var)[(index)]
#else
#   define GET_ITEM(var, index) (var).m128_f32[(index)]
#endif

template<typename T>
T Max(T a, T b)
{
    return (a > b) ? a : b;
}

template<typename T>
T Min(T a, T b)
{
    return (a < b) ? a : b;
}

typedef vector<float> TFloats;

static const float INF = 1e4f;

static void ConvertVector(const TFloats& floats, __m128** result, int* len)
{
    const int m = static_cast<int>(floats.size());
    *len = m / 4;
    if (m % 4)
    {
        ++(*len);
    }
    *result = (__m128*)_mm_malloc(sizeof(__m128)*(*len), 32);
    for (int i = 0; i < *len; ++i)
    {
        float values[] = { INF, INF, INF, INF };
        for (int j = 0; j < 4; ++j)
        {
            if (4 * i + j < m)
            {
                values[j] = floats[4 * i + j];
            }
        }
        (*result)[i] = _mm_load_ps(values);
    }
}

static float Median(const TFloats& floats)
{
    TFloats temp(floats);
    std::sort(temp.begin(), temp.end());
    return temp[temp.size() / 2];
}

struct KDTree
{
    int _len;
    __m128* _x;
    __m128* _y;
    bool _isLeaf;
    KDTree* _left;
    KDTree* _right;
    TIntVector _indices;
    float _minX;
    float _maxX;
    float _minY;
    float _maxY;

    KDTree(int depth, const TFloats& x, const TFloats& y, const TIntVector& indices)
        : _left(nullptr)
        , _right(nullptr)
    {
        if (depth == 9)
        {
            ConvertVector(x, &_x, &_len);
            ConvertVector(y, &_y, &_len);
            _indices = indices;
            _isLeaf = true;
            _minX = INF;
            _maxX = -INF;
            _minY = INF;
            _maxY = -INF;
            for (size_t i = 0; i < x.size(); ++i)
            {
                _minX = Min(_minX, x[i]);
                _maxX = Max(_maxX, x[i]);
                _minY = Min(_minY, y[i]);
                _maxY = Max(_maxY, y[i]);
            }
        }
        else
        {
            float separator = (depth & 1) ? Median(x) : Median(y);

            TFloats xLeft;
            TFloats yLeft;
            TIntVector indicesLeft;
            TFloats xRight;
            TFloats yRight;
            TIntVector indicesRight;

            for (size_t i = 0; i < x.size(); ++i)
            {
                bool toLeft = (depth & 1) ? (x[i] < separator) : (y[i] < separator);
                if (toLeft)
                {
                    xLeft.push_back(x[i]);
                    yLeft.push_back(y[i]);
                    indicesLeft.push_back(indices[i]);
                }
                else
                {
                    xRight.push_back(x[i]);
                    yRight.push_back(y[i]);
                    indicesRight.push_back(indices[i]);
                }
            }

            if (xLeft.size())
            {
                _left = new KDTree(depth + 1, xLeft, yLeft, indicesLeft);
            }
            if (xRight.size())
            {
                _right = new KDTree(depth + 1, xRight, yRight, indicesRight);
            }

            _isLeaf = false;
        }
    }

#ifndef _MSC_VER
    __attribute__((force_align_arg_pointer))
#endif
    void Solve(float x, float y, const __m128& x4, const __m128& y4, float* minDist, __m128* minMax, float* minMin, TIntVector* result) const
    {
        if (_isLeaf)
        {
            if (!_len)
            {
                return;
            }

            float d;
            if (x < _minX)
            {
                if (y < _minY)
                {
                    d = Min(_minX - x, _minY - y);
                }
                else if (y > _maxY)
                {
                    d = Min(_minX - x, y - _maxY);
                }
                else
                {
                    d = _minX - x;
                }
            }
            else if (x > _maxX)
            {
                if (y < _minY)
                {
                    d = Min(x - _maxX, _minY - y);
                }
                else if (y > _maxY)
                {
                    d = Min(x - _maxX, y - _maxY);
                }
                else
                {
                    d = x - _maxX;
                }
            }
            else
            {
                if (y < _minY)
                {
                    d = _minY - y;
                }
                else if (y > _maxY)
                {
                    d = y - _maxY;
                }
                else
                {
                    d = 0;
                }
            }

            if (d*d > *minDist)
            {
                return;
            }

            static const float EPS = 3e-7f;
            for (int j = 0; j < _len; ++j)
            {
                __m128 dx4 = _mm_sub_ps(_x[j], x4);
                dx4 = _mm_mul_ps(dx4, dx4);
                __m128 dy4 = _mm_sub_ps(_y[j], y4);
                dy4 = _mm_mul_ps(dy4, dy4);
                dx4 = _mm_add_ps(dx4, dy4);

                __m128 cmpMax = _mm_cmplt_ps(dx4, *minMax);
                for (int k = 0; k < 4; ++k)
                {
                    if (GET_ITEM(cmpMax, k))
                    {
                        if (GET_ITEM(dx4, k) < *minDist)
                        {
                            if (GET_ITEM(dx4, k) < *minMin)
                            {
                                result->clear();
                            }
                            *minDist = GET_ITEM(dx4, k);
                            *minMax = _mm_set1_ps(*minDist + EPS);
                            *minMin = *minDist - EPS;
                        }
                        result->push_back(_indices[4 * j + k]);
                    }
                }
            }
        }
        else
        {
            if (_left)
            {
                _left->Solve(x, y, x4, y4, minDist, minMax, minMin, result);
            }
            if (_right)
            {
                _right->Solve(x, y, x4, y4, minDist, minMax, minMin, result);
            }
        }
    }
};

void GenBig()
{
    FILE* fOut = fopen("big.txt", "w");
    static const int N = 100000;
    fprintf(fOut, "%d\n", N);
    for (int i = 0; i < N; ++i)
    {
        long double ldi = i;
        long double angle = ldi/N*2.0*M_PI;
        static const long double R = 10000.0;
        long double x = R*cos(angle);
        long double y = R*sin(angle);
        fprintf(fOut, "%Lf %Lf\n", x, y);
    }
    static const int M = 10000;
    fprintf(fOut, "%d\n", M);
    for (int i = 0; i < M; ++i)
    {
        long double ldi = i;
        long double angle = ldi/M*2.0*M_PI;
        static const long double R = 0.1;
        long double x = R*cos(angle);
        long double y = R*sin(angle);
        fprintf(fOut, "%Lf %Lf\n", x, y);
    }
    fclose(fOut);
}

int main()
{
#ifndef ONLINE_JUDGE
    GenBig();
    freopen("big.txt", "r", stdin);
    // freopen("input.txt", "r", stdin);
#endif

#ifndef __linux__
    SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
#endif

    int m;
    scanf("%d", &m);

    int size = m + 8;

    long double* xcd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    long double* ycd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    long double mx = 0.001;
    for (int i = 0; i < m; ++i)
    {
        scanf("%Lf%Lf", &xcd[i], &ycd[i]);
        mx = Max(mx, abs(xcd[i]));
        mx = Max(mx, abs(ycd[i]));
    }

    TIntVector indices(m);
    for (int i = 0; i < m; ++i)
    {
        indices[i] = i;
    }

    for (int i = 0; i < m; ++i)
    {
        swap(indices[i], indices[i + (rand() % (m - i))]);
    }

    TFloats xc(m);
    TFloats yc(m);
    for (int i = 0; i < m; ++i)
    {
        xc[i] = xcd[indices[i]] / mx;
        yc[i] = ycd[indices[i]] / mx;
    }

    KDTree* kdTree = new KDTree(0, xc, yc, indices);

    int n;
    scanf("%d", &n);
    vector<int> result;
    result.reserve(m);
    vector<int> result2;
    result2.reserve(m);
    for (int i = 0; i < n; ++i)
    {
        long double xd, yd;
        scanf("%Lf%Lf", &xd, &yd);

        float x = xd/mx;
        float y = yd/mx;
        __m128 x4 = _mm_set1_ps(x);
        __m128 y4 = _mm_set1_ps(y);

        result.clear();
        float minDist = INF*INF;
        __m128 minMax = _mm_set1_ps(minDist);
        float minMin = minDist;
        kdTree->Solve(x, y, x4, y4, &minDist, &minMax, &minMin, &result);

        if (result.size() < 1000)
        {
            result2.clear();
            long double mind = 1e15;
            long double mindMin = mind;
            long double mindMax = mind;
            for (auto realindex : result)
            {
                long double dist = Sqr(xd - xcd[realindex]) + Sqr(yd - ycd[realindex]);
                static const long double LDEPS = 1e-10;
                if (dist < mindMin)
                {
                    mind = dist;
                    mindMin = dist - LDEPS;
                    mindMax = dist + LDEPS;
                    result2.clear();
                }
                if (dist <= mindMax)
                {
                    result2.push_back(realindex);
                }
            }
        }
        else
        {
            result2 = result;
        }

        std::sort(result2.begin(), result2.end());

        Output(result2);
    }

    return 0;
}
