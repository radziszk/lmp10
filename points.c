#include "points.h"
#include <stdlib.h>

/*
	Powiększanie struktury punktów
*/

static int
realloc_pts_failed (points_t * pts, int size)
{
	double *iks =pts->x;
	double *igrek= pts->y; 
  int stan= (pts->x=realloc (pts->x, size * sizeof *pts->x)) == NULL
    || (pts->y=realloc (pts->y, size * sizeof *pts->y)) == NULL;
	if(stan)
	{
		free(pts->x);
		free(pts->y);
		free(igrek);
		free(iks);
	}
return stan;
}

/*
	Wczytanie struktury punktów z pliku
*/

int
read_pts_failed (FILE * inf, points_t * pts)
{
  int size;
  double x, y;

  if (pts->n == 0) {
    pts->x = malloc (100 * sizeof *pts->x);
    if (pts->x == NULL)
      return 1;
    pts->y = malloc (100 * sizeof *pts->y);
    if (pts->y == NULL) {
      free (pts->x);
      return 1;
    }
    size = 100;
  }
  else
    size = pts->n;

  while (fscanf (inf, "%lf %lf", &x, &y) == 2) {
    pts->x[pts->n] = x;
    pts->y[pts->n] = y;
    pts->n++;
    if (!feof (inf) && pts->n == size) {
      if (realloc_pts_failed (pts, 2 * size))
        return 1;
      else
        size *= 2;
    }
  }

  if (pts->n != size)
    if (realloc_pts_failed (pts, pts->n))
      return 1;

  return 0;
}
