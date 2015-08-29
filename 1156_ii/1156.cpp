#include <iostream>
#include <fstream>

#include <vector>
#include <stack>

using namespace std;

typedef vector<bool> TBoolVector;
typedef vector<TBoolVector> TBoolMatrix;

typedef vector<int> TIntVector;
typedef stack<int> TIntStack;

class TComponent
{
public:
    TIntVector _white;
    TIntVector _black;
    bool _first;
};

typedef vector<TComponent> TComponents;

void Impossible()
{
    cout << "IMPOSSIBLE" << endl;
}

bool Search(int leftWhite, int leftBlack, int index, TComponents& components)
{
    if (leftWhite < 0)
    {
        return false;
    }

    if (leftBlack < 0)
    {
        return false;
    }

    if (!leftWhite && !leftBlack)
    {
        return true;
    }

    components[index]._first = true;
    if (Search(leftWhite - components[index]._white.size(), leftBlack - components[index]._black.size(), index + 1, components))
    {
        return true;
    }

    components[index]._first = false;
    if (Search(leftWhite - components[index]._black.size(), leftBlack - components[index]._white.size(), index + 1, components))
    {
        return true;
    }

    return false;
}

int main()
{
#ifndef ONLINE_JUDGE
    ifstream fIn("input.txt");
    cin.rdbuf( fIn.rdbuf() );
#endif

    int n;
    int m;
    cin >> n >> m;

    const int n2 = n + n;

    TBoolMatrix matrix(n2, TBoolVector(n2, false));
    for (int i = 0; i < m; ++i)
    {
        int a;
        int b;
        cin >> a >> b;
        --a;
        --b;
        matrix[a][b] = true;
        matrix[b][a] = true;
    }

    TComponents components;
    TIntVector color(n2);
    for (int i = 0; i < n2; ++i)
    {
        if (0 == color[i])
        {
            components.push_back( TComponent() );
            TComponent& c = components.back();
            color[i] = 1;
            c._white.push_back(i);
            TIntStack s;
            s.push(i);
            while (!s.empty())
            {
                int top = s.top();
                s.pop();
                for (int j = 0; j < n2; ++j)
                {
                    if (matrix[top][j])
                    {
                        if (!color[j])
                        {
                            color[j] = -color[top];
                            if (1 == color[j])
                            {
                                c._white.push_back(j);
                            }
                            else
                            {
                                c._black.push_back(j);
                            }
                            s.push(j);
                        }
                    }
                }
            }
        }
    }

    bool failed = false;
    for (int i = 0; i < n2; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            if (matrix[i][j] && color[i] == color[j])
            {
                failed = true;
            }
        }
    }

    if (!failed)
    {
        if (Search(n, n, 0, components))
        {
            for (int i = 0; i < components.size(); ++i)
            {
                const TIntVector& v = (components[i]._first) ? components[i]._white : components[i]._black;
                for (int j = 0; j < v.size(); ++j)
                {
                    cout << v[j] + 1 << " ";
                }
            }
            cout << endl;
            for (int i = 0; i < components.size(); ++i)
            {
                const TIntVector& v = (components[i]._first) ? components[i]._black : components[i]._white;
                for (int j = 0; j < v.size(); ++j)
                {
                    cout << v[j] + 1 << " ";
                }
            }
            cout << endl;
        }
        else
        {
            Impossible();
        }
    }
    else
    {
        Impossible();
    }

    return 0;
}
