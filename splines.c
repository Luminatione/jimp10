#include "splines.h"

#include <stdlib.h>

#define MALLOC_FAILED( P, SIZE ) (((P)=malloc( (SIZE)*sizeof( *(P))))==NULL)

void freeSpline(spline_t* spline)
{
	if (spline->x == NULL)
		free(spline->x);
	if (spline->f == NULL)
		free(spline->f);
	if (spline->f1 == NULL)
		free(spline->f1);
	if (spline->f2 == NULL)
		free(spline->f2);
	if (spline->f3 == NULL)
		free(spline->f3);
}

int
alloc_spl(spline_t* spl, int n)
{
	spl->n = n;
	return MALLOC_FAILED(spl->x, spl->n)
		|| MALLOC_FAILED(spl->f, spl->n)
		|| MALLOC_FAILED(spl->f1, spl->n)
		|| MALLOC_FAILED(spl->f2, spl->n)
		|| MALLOC_FAILED(spl->f3, spl->n);
}

int
read_spl(FILE* inputFile, spline_t* spl)
{
	int i;
	if (fscanf(inputFile, "%d", &(spl->n)) != 1 || spl->n < 0)
		return 1;

	if (alloc_spl(spl, spl->n))
		return 1;

	for (i = 0; i < spl->n; i++)
		if (fscanf
		(inputFile, "%lf %lf %lf %lf %lf", spl->x + i, spl->f + i, spl->f1 + i,
			spl->f2 + i, spl->f3 + i) != 5)
			return 1;

	return 0;
}

void
write_spl(spline_t* spl, FILE* outputFile)
{
	int i;
	fprintf(outputFile, "%d\n", spl->n);
	for (i = 0; i < spl->n; i++)
		fprintf(outputFile, "%g %g %g %g %g\n", spl->x[i], spl->f[i], spl->f1[i],
			spl->f2[i], spl->f3[i]);
}

double
value_spl(spline_t* spl, double x)
{
	int i;
	double dx;

	for (i = spl->n - 1; i > 0; i--)
		if (spl->x[i] < x)
			break;

	dx = x - spl->x[i];

	return spl->f[i]
		+ dx * spl->f1[i]
		+ dx * dx / 2 * spl->f2[i]
		+ dx * dx * dx / 6 * spl->f3[i];
}
