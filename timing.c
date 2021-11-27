#include <time.h>

void get_wtime(double *wt)
{
  clock_t t = clock();
  *wt = (double)t /  CLOCKS_PER_SEC;
}

void get_wtime_(double *wt)
{
	get_wtime(wt);
}
