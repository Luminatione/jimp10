#ifndef POINTS_H
#define POINTS_H

#include <stdio.h>

typedef struct {
		int n;
		double *x;
		double *y;
} points_t;

void freePoints(points_t* points);

int read_pts_failed ( FILE* inf, points_t *pts);

#endif
