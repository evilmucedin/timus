#include <cstdio>

#include <vector>

using namespace std;

typedef vector<int> TIntVector;

struct TVertex
{
    int _information;
    TIntVector _edges;
};

typedef vector<TVertex> TVertexes;

int main()
{
    int n;
    scanf("%d", &n);

    TVertexes vertexes(n);
    for (int i = 0; i < n; ++i)
    {
        scanf("%d", &vertexes[i]._information);
    }
    for (int i = 1; i < n; ++i)
    {
        int a;
        int b;
        scanf("%d%d", &a, &b);
        --a;
        --b;
        vertexes[a]._edges.push_back(b);
        vertexes[b]._edges.push_back(a);
    }

    int m;
    scanf("%d", &m);
    for (int i = 0; i < m; ++i)
    {
        int t;
        int v;
        scanf("%d%d", &t, &v);
        --v;
        if (1 == t)
        {
            for (const auto& e : vertexes[v]._edges)
            {
                vertexes[e]._information = (vertexes[e]._information + vertexes[v]._information) % 1000000007;
            }
        }
        else
        {
            printf("%d\n", vertexes[v]._information);
        }
    }

    return 0;
}
