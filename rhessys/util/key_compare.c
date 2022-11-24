#include <stdlib.h>
#include "rhessys.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
/*------------------------------
sort heights in descending order
-------------------------------*/  
int key_compare( void * e1,  void *e2 )
{
	/*------------------------------------------------------*/
	/*	Local Function Definition. 							*/
	/*------------------------------------------------------*/
	double v1,v2;	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
		v1 = ((struct layer_object *)e1)->height;
		v2 = ((struct layer_object *)e2)->height;

	return (v1>v2) ? -1 : (v1<v2) ? 1: 0;
}/*end key_compare.c*/
//05202022LML compare double
//int close_enough(double a, double b)
//{
//    if (fabs(a - b) <= DBL_EPSILON * fmax(fabs(a), fabs(b))) return 1;
//    return 0;
//}
