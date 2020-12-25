#include "main.h"


int main(){
	int i;
	tlist list;
	point p[5];
	tlistinit(&list);
	p[0].x = 0.0; p[0].y=0.0;
	p[1].x = 2.0; p[1].y=0.0;
	p[2].x = 0.0; p[2].y=2.0;
	p[3].x = 2.0; p[3].y=2.0;
	p[4].x = 1.0; p[4].y=1.0;
	BowyerWatson( &list , p , 5 );
	for(i=0; i<list.total; ++i){
		printf("%d:(%f %f), (%f %f), (%f %f)\n", i, list.item[i].a.x, list.item[i].a.y\
							  , list.item[i].b.x, list.item[i].b.y\
							  , list.item[i].c.x, list.item[i].c.y);
	}
	tlistfree(&list);
}

void BowyerWatson( tlist *t , point *p , int np ){
	size_t i , j , k , l, n,m;
	bool flag;
	ilist bad_triangle_id;
	elist polygon;
	edge e[3];
	point P3[3];
	triangle temp;
	//super triangle
	P3[0].x = -500;
	P3[0].y = -200;
	P3[1].x = 500;
	P3[1].y = -200;
	P3[2].x = 0;
	P3[2].y = 300;
	new_triangle_p( &temp , P3[0] , P3[1] , P3[2], &flag );
	tlistpush(t , temp );

	ilistinit(&bad_triangle_id);
	elistinit(&polygon);
	for( l = 0 ; l < np ; ++l){
		for( i = 0 , n = t->total ; i < n ; ++i ){
			if( circum_contain( t->item[i] , p[l] )){
				ilistpush( &bad_triangle_id , i );
			}
		}
		for( i = 0 , n = bad_triangle_id.total ; i < n ; ++i ){
			to_edge(t->item[bad_triangle_id.item[i]] , e );
			for( j = 0; j < 3 ; ++j ){
				for ( k = 0; k < n; ++k ){
					if( i != k && edge_contain(t->item[bad_triangle_id.item[k]],e[j]) ){
						break;
					}
				}
				if( k == n ){
					elistpush( &polygon, e[j] );
				}
			}
		}
		k=0;
		for( i = 0 , n = bad_triangle_id.total ; i < n; ++i){
			new_triangle_p( &temp, p[l] , polygon.item[i].a , polygon.item[i].b , &flag);
			if(flag){
				tlistset(t, temp, bad_triangle_id.item[i]-k);
			}else{
				tlistdelete( t, bad_triangle_id.item[i]-k);
				k++;
			}
		}
		for( n = polygon.total ; i < n; ++i){
			new_triangle_p( &temp, p[l] , polygon.item[i].a , polygon.item[i].b, &flag);
			if(flag)tlistpush(t, temp);
		}
		polygon.total=0;
		bad_triangle_id.total=0;
	}
	for(j = 0 ; j<3 ; j++){
		for(i = t->total; i>0 ; --i){
			if(point_contain(t->item[i], P3[j])){
				tlistdelete(t,i);
			}		
		}
		if(point_contain(t->item[0], P3[j])){
			t->total--;
			for(i=0; i<t->total; i++){
				new_triangle_t(&(t->item[i]), t->item[i+1]);
			}
		}	
	}
	elistfree( &polygon );
	ilistfree( &bad_triangle_id );
}

void check(triangle *a, bool *flag){
	float d = (a->b.x - a->a.x) * (a->c.y - a->a.y)\
		- (a->b.y - a->a.y) * (a->c.x - a->a.x);
	if(d==0.0){
		printf("three point on the same line!\n");
		*flag = 0;
		return;
	}else if(d < 0){
		float temp;
		temp = a->a.x;
		a->a.x = a->b.x;
		a->b.x = temp;
		temp = a->a.y;
		a->a.y = a->b.y;
		a->b.y = temp;
	}
	*flag=1;
}

void circum(triangle *a){
	float D;
	D = 2 * ( a->a.x * ( a->b.y - a->c.y )\
	        + a->b.x * ( a->c.y - a->a.y )\
	        + a->c.x * ( a->a.y - a->b.y ) );
	a->center.x = ( ( a->a.x * a->a.x + a->a.y * a->a.y )*( a->b.y - a->c.y )\
		      + ( a->b.x * a->b.x + a->b.y * a->b.y )*( a->c.y - a->a.y )\
		      + ( a->c.x * a->c.x + a->c.y * a->c.y )*( a->a.y - a->b.y ) ) / D;
	a->center.y = ( ( a->a.x * a->a.x + a->a.y * a->a.y )*( a->c.x - a->b.x )\
		      + ( a->b.x * a->b.x + a->b.y * a->b.y )*( a->a.x - a->c.x )\
		      + ( a->c.x * a->c.x + a->c.y * a->c.y )*( a->b.x - a->a.x ) ) / D;
	a->r_sqr = (a->center.x - a->a.x)*(a->center.x-a->a.x)+(a->center.y-a->a.y)*(a->center.y-a->a.y);
}

void to_edge( triangle a , edge*l ){
	l[0].a.x = a.a.x;
	l[0].a.y = a.a.y;
	l[0].b.x = a.b.x;
	l[0].b.y = a.b.y;
	l[1].a.x = a.b.x;
	l[1].a.y = a.b.y;
	l[1].b.x = a.c.x;
	l[1].b.y = a.c.y;
	l[2].a.x = a.c.x;
	l[2].a.y = a.c.y;
	l[2].b.x = a.a.x;
	l[2].b.y = a.a.y;
}

void  new_triangle_p(triangle *a, point p1, point p2, point p3, bool *flag){
	a->a.x = p1.x;
	a->a.y = p1.y;
	a->b.x = p2.x;
	a->b.y = p2.y;
	a->c.x = p3.x;
	a->c.y = p3.y;
	check(a, flag);
	circum(a);
}

void  new_triangle_t(triangle *a, triangle b){
	a->a.x = b.a.x;
	a->a.y = b.a.y;
	a->b.x = b.b.x;
	a->b.y = b.b.y;
	a->c.x = b.c.x;
	a->c.y = b.c.y;
	a->center.x = b.center.x;
	a->center.y = b.center.y;
	a->r_sqr = b.r_sqr;
}

void new_edge(edge *a, point p1, point p2){
	a->a.x = p1.x;
	a->a.y = p1.y;
	a->b.x = p2.x;
	a->b.y = p2.y;	
}	

bool circum_contain(triangle a , point x){
	return ( ( a.center.x - x.x ) * ( a.center.x - x.x )\
			+ ( a.center.y - x.y ) * ( a.center.y - x.y ) ) < a.r_sqr;
}

bool edge_contain(triangle a, edge l){
	return ( ( a.a.x == l.a.x ) && ( a.a.y == l.a.y ) && ( a.b.x == l.b.x ) && ( a.b.y == l.b.y ) )\
	    || ( ( a.b.x == l.a.x ) && ( a.b.y == l.a.y ) && ( a.a.x == l.b.x ) && ( a.a.y == l.b.y ) )\
	    || ( ( a.b.x == l.a.x ) && ( a.b.y == l.a.y ) && ( a.c.x == l.b.x ) && ( a.c.y == l.b.y ) )\
	    || ( ( a.c.x == l.a.x ) && ( a.c.y == l.a.y ) && ( a.b.x == l.b.x ) && ( a.b.y == l.b.y ) )\
	    || ( ( a.c.x == l.a.x ) && ( a.c.y == l.a.y ) && ( a.a.x == l.b.x ) && ( a.a.y == l.b.y ) )\
	    || ( ( a.a.x == l.a.x ) && ( a.a.y == l.a.y ) && ( a.c.x == l.b.x ) && ( a.c.y == l.b.y ) );
}

bool point_contain(triangle a, point p ){
	return ( ( a.a.x == p.x ) && ( a.a.y == p.y ) )\
	    || ( ( a.b.x == p.x ) && ( a.b.y == p.y ) )\
	    || ( ( a.c.x == p.x ) && ( a.c.y == p.y ) );
}

bool nextto(triangle a, triangle b){
	return ( ( a.b.x == b.a.x ) && ( a.b.y == b.a.y ) && ( a.a.x == b.b.x ) || ( a.a.y == b.b.y ) )\
	    || ( ( a.a.x == b.a.x ) && ( a.a.y == b.a.y ) && ( a.c.x == b.b.x ) && ( a.c.y == b.b.y ) )\
	    || ( ( a.c.x == b.a.x ) && ( a.c.y == b.a.y ) && ( a.b.x == b.b.x ) && ( a.b.y == b.b.y ) )\
	    || ( ( a.b.x == b.b.x ) && ( a.b.y == b.b.y ) && ( a.a.x == b.c.x ) && ( a.a.y == b.c.y ) )\
	    || ( ( a.a.x == b.b.x ) && ( a.a.y == b.b.y ) && ( a.c.x == b.c.x ) && ( a.c.y == b.c.y ) )\
	    || ( ( a.c.x == b.b.x ) && ( a.c.y == b.b.y ) && ( a.b.x == b.c.x ) && ( a.b.y == b.c.y ) )\
	    || ( ( a.b.x == b.c.x ) && ( a.b.y == b.c.y ) && ( a.a.x == b.a.x ) && ( a.a.y == b.a.y ) )\
	    || ( ( a.a.x == b.c.x ) && ( a.a.y == b.c.y ) && ( a.c.x == b.a.x ) && ( a.c.y == b.a.y ) )\
	    || ( ( a.c.x == b.c.x ) && ( a.c.y == b.c.y ) && ( a.b.x == b.a.x ) && ( a.b.y == b.a.y ) );
}

void tlistinit(tlist *list){
	list->capacity = INITIAL;
	list->total = 0;
	list->item = (triangle*)malloc( sizeof(triangle) * list->capacity);
}

void elistinit(elist *list){
	list->capacity = INITIAL;
	list->total = 0;
	list->item = (edge*)malloc( sizeof(edge) * list->capacity);
}

void ilistinit(ilist *list){
	list->capacity = INITIAL;
	list->total = 0;
	list->item = (int*)malloc( sizeof(int) * list->capacity);
}

void tlistpush(tlist *list, triangle a){
	if(list){
		if( list->capacity == list->total){
			list->capacity += INITIAL;
			list->item = (triangle*)realloc( list->item, sizeof(triangle) * list->capacity );
		}
		new_triangle_t( &( list->item[ list->total++ ] ) , a);
	}
}

void elistpush(elist *list, edge a){
	if(list){
		if( list->capacity == list->total){
			list->capacity += INITIAL;
			list->item = (edge*)realloc( list->item, sizeof(edge) * list->capacity );
		}
		new_edge( &( list->item[ list->total++ ] ) , a.a , a.b );
	}
}

void ilistpush(ilist *list, int a){
	if(list){
		if( list->capacity == list->total){
			list->capacity += INITIAL;
			list->item = ( int* )realloc( list->item, sizeof(int) * list->capacity );
		}
		list->item[ list->total++ ] = a;
	}
}

void tlistset(tlist *list, triangle a, int index){
	if( (index >= 0) && (index < list->total ) ){
		new_triangle_t( &( list->item[ index ] ) , a );
	}
}
	
void elistset(elist *list, edge a, int index){
	if( (index >= 0) && (index < list->total ) ){
		new_edge( &( list->item[ index ] ) , a.a , a.b );
	}
}
void ilistset(ilist *list, int a, int index){
	if( (index >= 0) && (index < list->total ) ){
		list->item[ index ] = a;
	}
}

void tlistdelete(tlist *list , int index){
	int i;
	if( (index >= 0) && (index < list->total ) ){
		list->total--;
		for(i = index; i < list->total; ++i){
			new_triangle_t(&(list->item[i]), list->item[i+1]);
		}
	}
}

void elistdelete(elist *list , int index){
	int i;
	if( (index >= 0) && (index < list->total ) ){
		list->total--;
		for(i = index; i < list->total; ++i){
			new_edge(&(list->item[i]), list->item[i+1].a , list->item[i+1].b );
		}
	}
}

void ilistdelete(ilist *list , int index){
	int i;
	if( (index >= 0) && (index < list->total ) ){
		list->total--;
		for(i = index; i < list->total; ++i){
			list->item[i] = list->item[i+1];
		}	
	}
}

void tlistfree(tlist *list){
	free( list->item );
}

void elistfree(elist *list){
	free( list->item );
}

void ilistfree(ilist *list){
	free( list->item );
}
