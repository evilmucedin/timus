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

#define _CRT_SECURE_NO_WARNINGS

#pragma GCC target("sse4.2")
#pragma GCC optimize("O3")

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

__m128 x4;
__m128 y4;
__m128* xc4;
__m128* yc4;
int* indexes;

struct TThreadParam
{
    int begin;
    int end;
    TIntVector result;
    HANDLE thread;
};

 __attribute__((force_align_arg_pointer))
DWORD WINAPI Nop(LPVOID param)
{
    TThreadParam* tParam = (TThreadParam*)param;

    static const float INF = 1e16f;
    float min = INF;
    __m128 minMax = _mm_set1_ps(INF);
    float minMin = INF;

    tParam->result.clear();
    for (int j = tParam->begin; j < tParam->end; ++j)
    {
        __m128 dx4 = _mm_sub_ps(xc4[j], x4);
        dx4 = _mm_mul_ps(dx4, dx4);
        __m128 dy4 = _mm_sub_ps(yc4[j], y4);
        dy4 = _mm_mul_ps(dy4, dy4);
        dx4 = _mm_add_ps(dx4, dy4);

        __m128 cmpMax = _mm_cmplt_ps(dx4, minMax);
        for (int k = 0; k < 4; ++k)
        {
            if (cmpMax[k])
            {
                if (dx4[k] < min)
                {
                    static const float EPS = 3e-7;
                    if (dx4[k] < minMin)
                    {
                        tParam->result.clear();
                    }
                    min = dx4[k];
                    minMax = _mm_set1_ps(min + EPS);
                    minMin = min - EPS;
                }
                tParam->result.push_back(4 * j + k);
            }
        }
    }

    return 0;
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
    long double mx = 0.001;
    for (int i = 0; i < m; ++i)
    {
        scanf("%Lf%Lf", &xcd[i], &ycd[i]);
        mx = std::max(mx, abs(xcd[i]));
        mx = std::max(mx, abs(ycd[i]));
    }

    indexes = new int[size];
    for (int i = 0; i < size; ++i)
    {
        indexes[i] = i;
    }

    for (int i = 0; i < m; ++i)
    {
        swap(indexes[i], indexes[i + (rand() % (m - i))]);
    }

    float* xc = (float*)_mm_malloc(sizeof(float)*size, 32);
    float* yc = (float*)_mm_malloc(sizeof(float)*size, 32);
    for (int i = 0; i < size; ++i)
    {
        xc[i] = 1e5f;
        yc[i] = 1e5f;
    }
    for (int i = 0; i < m; ++i)
    {
        xc[i] = xcd[indexes[i]]/mx;
        yc[i] = ycd[indexes[i]]/mx;
    }

    const int mLen = m / 4 + 1;
    xc4 = (__m128*)_mm_malloc(sizeof(__m128)*(mLen), 32);
    yc4 = (__m128*)_mm_malloc(sizeof(__m128)*(mLen), 32);
    for (int i = 0; i < mLen; ++i)
    {
        xc4[i] = _mm_load_ps(xc + 4*i);
        yc4[i] = _mm_load_ps(yc + 4*i);
    }

    static const size_t NTHREADS = 1;
    vector<TThreadParam> params(NTHREADS);
    int begin = 0;
    int step = mLen/NTHREADS;
    for (int i = 0; i < NTHREADS; ++i)
    {
        params[i].begin = begin;
        params[i].end = begin + step;
        begin += step;
    }
    params[NTHREADS - 1].end = mLen;

    int n;
    scanf("%d", &n);
    vector<int> result2;
    result2.reserve(m);
    for (int i = 0; i < n; ++i)
    {
        long double xd, yd;
        scanf("%Lf%Lf", &xd, &yd);

        float x, y;
        x = xd/mx;
        y = yd/mx;
        x4 = _mm_set1_ps(x);
        y4 = _mm_set1_ps(y);

#if !defined(__linux__) && defined(MT_1369)
	/*	
	The Timus machine has four cores.	
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	if (4 == si.dwNumberOfProcessors)
	{
            return 0;	
	}
	*/

        for (int j = 0; j < NTHREADS; ++j)
        {
            params[j].thread = CreateThread(0, 0, Nop, &(params[j]), 0, 0);
            if (!params[j].thread)
            {
                return 0;
            }
            SetThreadPriority(params[j].thread, THREAD_PRIORITY_TIME_CRITICAL);
            SetThreadAffinityMask(params[j].thread, 1 << j);
        }

        for (int j = 0; j < NTHREADS; ++j)
        {
            WaitForSingleObject(params[j].thread, INFINITE);
        }
#else
        for (size_t iThread = 0; iThread < NTHREADS; ++iThread)
        {
            Nop(&(params[iThread]));
        }
#endif

        result2.clear();
        long double mind = 1e15;
        long double mindMin = mind;
        long double mindMax = mind;
        for (int iThread = 0; iThread < NTHREADS; ++iThread)
        {
            for (auto index : params[iThread].result)
            {
                int realindex = indexes[index];
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
        std::sort(result2.begin(), result2.end());

        Output(result2);
    }

    return 0;
}
