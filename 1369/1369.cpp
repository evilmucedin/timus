#define _CRT_SECURE_NO_WARNINGS
#define M_PI       3.14159265358979323846

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

union TVector
{
    __m128d v;
    double d[2];
};

typedef vector<double> TDoubles;

static const double INF = 1e4f;

static void ConvertVector(const TDoubles& floats, TVector** result, int* len)
{
    const int m = static_cast<int>(floats.size());
    *len = m / 2;
    if (m % 2)
    {
        ++(*len);
    }
    *result = (TVector*)_mm_malloc(sizeof(TVector)*(*len), 32);
    for (int i = 0; i < *len; ++i)
    {
        double values[] = { INF, INF };
        for (int j = 0; j < 2; ++j)
        {
            if (2 * i + j < m)
            {
                values[j] = floats[2 * i + j];
            }
        }
        (*result)[i].v = _mm_load_pd(values);
    }
}

static double Median(const TDoubles& floats)
{
    TDoubles temp(floats);
    std::sort(temp.begin(), temp.end());
    return temp[temp.size() / 2];
}

struct KDTree
{
    int _len;
    TVector* _x;
    TVector* _y;
    bool _isLeaf;
    KDTree* _left;
    KDTree* _right;
    TIntVector _indices;
    double _minX;
    double _maxX;
    double _minY;
    double _maxY;

    KDTree(int depth, const TDoubles& x, const TDoubles& y, const TIntVector& indices)
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
            double separator = (depth & 1) ? Median(x) : Median(y);

            TDoubles xLeft;
            TDoubles yLeft;
            TIntVector indicesLeft;
            TDoubles xRight;
            TDoubles yRight;
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
    void Solve(double x, double y, const TVector& x2, const TVector& y2, double* minDist, TVector* minMax, double* minMin, TIntVector* result) const
    {
        if (_isLeaf)
        {
            if (!_len)
            {
                return;
            }

            double d;
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

            static const double EPS = 1e-13;
            for (int j = 0; j < _len; ++j)
            {
                TVector dx2;
                dx2.v = _mm_sub_pd(_x[j].v, x2.v);
                dx2.v = _mm_mul_pd(dx2.v, dx2.v);
                TVector dy2;
                dy2.v = _mm_sub_pd(_y[j].v, y2.v);
                dy2.v = _mm_mul_pd(dy2.v, dy2.v);
                dx2.v = _mm_add_pd(dx2.v, dy2.v);

                TVector cmpMax;
                cmpMax.v = _mm_cmplt_pd(dx2.v, minMax->v);
                for (int k = 0; k < 2; ++k)
                {
                    if (cmpMax.d[k])
                    {
                        if (dx2.d[k] < *minDist)
                        {       
                            if (dx2.d[k] < *minMin)
                            {
                                result->clear();
                            }
                            *minDist = dx2.d[k];
                            minMax->v = _mm_set1_pd(*minDist + EPS);
                            *minMin = *minDist - EPS;
                        }
                        result->push_back(_indices[2 * j + k]);
                    }
                }
            }
        }
        else
        {
            if (_left)
            {
                _left->Solve(x, y, x2, y2, minDist, minMax, minMin, result);
            }
            if (_right)
            {
                _right->Solve(x, y, x2, y2, minDist, minMax, minMin, result);
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
    // GenBig();
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
    long double mx = 1.0;
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

    TDoubles xc(m);
    TDoubles yc(m);
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

        double x = xd/mx;
        double y = yd/mx;
        TVector x2;
        x2.v = _mm_set1_pd(x);
        TVector y2;
        y2.v = _mm_set1_pd(y);

        result.clear();
        double minDist = 1e8;
        TVector minMax;
        minMax.v = _mm_set1_pd(minDist);
        double minMin = minDist;
        kdTree->Solve(x, y, x2, y2, &minDist, &minMax, &minMin, &result);

        if (result.size() < 10)
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
