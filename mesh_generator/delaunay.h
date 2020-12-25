#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#define INITIAL 10
typedef struct {
	float x;
	float y;
}point;

typedef struct {
	point a, b, c;
	point center;
	float r_sqr;
}triangle;

typedef struct {
	point a, b;
}edge;

typedef struct {
	int *item;
	int capacity;
	int total;
}ilist;

typedef struct triangle_list{
	triangle *item;
	int capacity;
	int total;
}tlist;

typedef struct edge_list{
	edge *item;
	int capacity;
	int total;
}elist;

void BowyerWatson( tlist *t , point *p , int np );

void check(triangle *a, bool *flag);
void circum(triangle *a);
void to_edge(triangle a, edge *l);
void new_triangle_p(triangle *a, point p1, point p2, point p3, bool *flag);
void new_triangle_t(triangle *a, triangle b);
void new_edge(triangle *a, point p1, point p2 );
bool circum_contain(triangle a, point x);
bool edge_contain(triangle a , edge l );
bool point_contain( triangle a, point x );
bool nextto(triangle a, triangle b );

void tlistinit(tlist *list);
void elistinit(elist *list);
void ilistinit(ilist *list);
void tlistpush(tlist *list, triangle a);
void elistpush(elist *list, edge a);
void ilistpush(ilist *list, int a);
void tlistset(tlist *list, triangle a, int index);
void elistset(elist *list, edge a, int index);
void ilistset(ilist *list, int a, int index);
void tlistdelete(tlist *list , int index);
void elistdelete(elist *list , int index);
void ilistdelete(ilist *list , int index);
void tlistfree(tlist *list);
void elistfree(elist *list);
void ilistfree(ilist *list);
