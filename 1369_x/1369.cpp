#define _CRT_SECURE_NO_WARNINGS

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

#ifndef _MSC_VER
#   define GET_ITEM2(var, index) (var)[(index)]
#else
#   define GET_ITEM2(var, index) (var).m128d_f64[(index)]
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

static const float INF = 1e6f;

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
    size_t n = temp.size() / 2;
    nth_element(temp.begin(), temp.begin() + n, temp.end());
    return temp[n];
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
#endif

#ifndef __linux__
    SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
#endif

    int m;
    scanf("%d", &m);

    int size = m + 8;

    long double* xcd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    long double* ycd = (long double*)_mm_malloc(sizeof(long double)*size, 32);
    __m128d* xyd = (__m128d*)_mm_malloc(sizeof(__m128d)*size, 32);
    for (int i = 0; i < m; ++i)
    {
        scanf("%Lf%Lf", &xcd[i], &ycd[i]);
        xyd[i] = _mm_set_pd(xcd[i], ycd[i]);
    }

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
        __m128d xy = _mm_set_pd(xd, yd);

        long double mind = 1e100;
        long double mindMin = mind;
        long double mindMax = mind;
        for (int j = 0; j < m; ++j)
        {
            __m128d xyd2 = _mm_sub_pd(xyd[j], xy);
            xyd2 = _mm_mul_pd(xyd2, xyd2);
            long double dist = GET_ITEM2(xyd2, 0) + GET_ITEM2(xyd2, 1);
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
                result2.push_back(j);
            }
        }

        Output(result2);
    }

    return 0;
}