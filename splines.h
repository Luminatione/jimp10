#ifndef SPLINES_H
#define SPLINES_H

#include <stdio.h>

typedef struct {
		int n;
		double *x;
		double *f;
		double *f1;
		double *f2;
		double *f3;
} spline_t;

void freeSpline(spline_t* spline);

int alloc_spl( spline_t *spl, int n );

int  read_spl ( FILE *inputFile,  spline_t *spl );

void  write_spl ( spline_t *spl, FILE * outputFile );

double value_spl( spline_t *spl, double x);

#endif
