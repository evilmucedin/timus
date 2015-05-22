#include	<stdio.h>

#define	MAXE	1011
#define	MAXF	1011

int cache[MAXE][MAXF];

inline int Min(int a, int b)
{
	return (a < b) ? a : b;
}

inline int Max(int a, int b)
{
	return (a > b) ? a : b;
}

int Solve(int eggs, int floor)
{
	int e, f;
	for (f = 0; f <= floor; ++f)
		cache[0][f] = 1000000;
	for (f = 1; f <= floor; ++f)
		cache[1][f] = f;
	for (e = 1; e <= eggs; ++e)
		cache[e][1] = 1;
	for (e = 1; e <= eggs; ++e)
		cache[e][0] = 0;

	for (e = 2; e <= eggs; ++e)
		for (f = 2; f <= floor; ++f)
		{
			int min = 1000000;
			for (int k = 1; k <= f; ++k)
			{
				int now = Max(cache[e - 1][k - 1], cache[e][f - k]) + 1;
				min = Min(min, now);
			}
			cache[e][f] = min;
		}

	return cache[eggs][floor];
}

int main(void)
{
	Solve(1003, 1003);
	int eggs, floor;
	while (scanf("%d %d", &eggs, &floor) == 2)
	{
		if (eggs == 0)
			break;
		if (eggs == 1)
			printf("%d\n", floor);
		else
			printf("%d\n", cache[eggs][floor]);
	}

	return 0;
}
