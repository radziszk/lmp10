#ifndef POINTS_H
#define POINTS_H

#include <stdio.h>

/*
	Struktura przechowująca punkty wczytane z pliku
*/

typedef struct {
		int n;
		double *x;
		double *y;
} points_t;


int read_pts_failed ( FILE* inf, points_t *pts);

#endif
