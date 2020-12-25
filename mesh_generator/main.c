#include "delaunay.h"
#include <math.h>
#include <time.h>

// one for yes, zero for no
#define UNIFORM 1
#define D 0.08

// model for fem
// _________
// |       |
// |       |
// |       |
// |       |__________
// |                 |
// |                 |
// |                 |
// |_________________|
//

int main(){
	int i , np = 100, j;
	tlist list;
	point *p;
	bool flag;
	int *ind;
	FILE *fp1, *fp2;
	float d;
	p = ( point* )malloc( np * sizeof(point) );
	ind = ( int* )malloc( 3*np*sizeof(point) );
	srand(time(NULL));
	p[0].x=0.0;
	p[0].y=0.0;
	p[1].x=1.0;
	p[1].y=0.0;
	p[2].x=0.0;
	p[2].y=1.0;
	p[3].x=1.0;
	p[3].y=1.0;
	p[4].x=0.5;
	p[4].y=0.0;
	p[5].x=0.0;
	p[5].y=0.5;
	p[6].x=1.0;
	p[6].y=0.5;
	p[7].x=0.5;
	p[7].y=1.0;
	p[8].x=0.25;
	p[8].y=0.0;
	p[9].x=0.0;
	p[9].y=0.25;
	p[10].x=0.75;
	p[10].y=0.0;
	p[11].x=0.75;
	p[11].y=0.0;
	p[12].x=1.0;
	p[12].y=0.25;
	p[13].x=.25;
	p[13].y=1.0;
	p[14].x=1.0;
	p[14].y=0.75;
	p[15].x=0.75;
	p[15].y=1.0;

	for(i=15; i<np; i++){
		if(UNIFORM){
			flag=1;
			while(flag){
				p[i].x=(float)rand()/RAND_MAX;
				p[i].y=(float)rand()/RAND_MAX;
				for(j=0; j<i; j++){
					d=(p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y);
					if(d<D*D){
						break;
					}
				}
				if(j==i){
					flag=0;
				}
			}	

		}else{	
			p[i].x=(float)rand()/RAND_MAX;
			p[i].y=(float)rand()/RAND_MAX;
		}
	}
	for(i=0; i<np; i++){
		printf("%f %f\n",p[i].x, p[i].y);
	}
	tlistinit(&list);
	BowyerWatson( &list , p , np );
	printf("total=%d\n",list.total);
	fp1 = fopen("triangle_x.txt", "w");
	fp2 = fopen("triangle_y.txt", "w");
	for(i=0; i<list.total; i++){
		fprintf(fp1, "%f ",list.item[i].a.x);	
		fprintf(fp1, "%f ",list.item[i].b.x);	
		fprintf(fp1, "%f\n",list.item[i].c.x);;
		fprintf(fp2, "%f ",list.item[i].a.y);	
		fprintf(fp2, "%f ",list.item[i].b.y);	
		fprintf(fp2, "%f\n",list.item[i].c.y);	
	}
	fclose(fp1);
	fclose(fp2);
	free(ind);
	free(p);
	tlistfree(&list);
}

