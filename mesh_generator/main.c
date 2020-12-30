#include "delaunay.h"
#include <math.h>
#include <time.h>

// one for yes, zero for no
#define UNIFORM 1
#define Np 150
#define D 1.7*(0.75/M_PI/Np)
//#define Nb ((int)(0.5/D)*4)
#define Nb 12
#define dx (1.0/Nb)

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
	int i , np = Np, j, k;
	tlist list;
	point *p;
	bool flag;
	FILE *fp1, *fp2, *fp3;
	float d, cx, cy;
	int count;
	p = ( point* )malloc( np * sizeof(point) );
	srand(time(NULL));
	p[0].x=0.0;   p[0].y=0.0;
	p[1].x=1.0;   p[1].y=0.0;
	p[2].x=0.0;   p[2].y=1.0;
	p[3].x=0.5;   p[3].y=0.5;
	p[4].x=0.5;   p[4].y=1.0;
	p[5].x=1.0;   p[5].y=0.5;
	i=6;

	for(j=1; j<Nb; j++){
		p[i].x = dx*j; p[i].y=0.0;
		i++;
		p[i].x = 0.0;    p[i].y=dx*j;
		i++;
	}
	for(j=1; j<0.5*Nb; j++){
		p[i].x = dx*j; p[i].y=1.0;
		i++;
		p[i].x = 1.0;      p[i].y=dx*j;
		i++;
	}
	for(j=1; j<0.5*Nb; j++){
		p[i].x = dx*j+0.5; p[i].y=0.5;
		i++;
		p[i].x = 0.5;        p[i].y=dx*j+0.5;
		i++;
	}

	for(; i<np; i++){
		if(UNIFORM){
			flag=1;
			while(flag){
				p[i].x=(float)rand()/RAND_MAX;
				p[i].y=(float)rand()/RAND_MAX;
				for(j=0; j<i; j++){
					if(p[i].x>=0.5&&p[i].y>=0.5){
						break;
					}
					d=(p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y);
					if(d<D){
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
	fp1 = fopen("point.txt","w");
	for(i=0; i<np; i++){
		fprintf(fp1, "%f %f\n",p[i].x, p[i].y);
	}
	fclose(fp1);

	tlistinit(&list);
	BowyerWatson( &list , p , np );
	for(i=list.total; i>=0;--i){
		cx=(list.item[i].a.x+list.item[i].b.x+list.item[i].c.x)/3;	
		cy=(list.item[i].a.y+list.item[i].b.y+list.item[i].c.y)/3;
		if(cx>0.5&&cy>0.5){
			tlistdelete( &list, i );
		}
	}
	printf("total=%d\n",list.total);
	fp1 = fopen("triangle_x.txt", "w");
	fp2 = fopen("triangle_y.txt", "w");
	fp3 = fopen("tri_index.txt", "w");
	for(i=0; i<list.total; i++){
		fprintf(fp1, "%f ",list.item[i].a.x);	
		fprintf(fp1, "%f ",list.item[i].b.x);	
		fprintf(fp1, "%f\n",list.item[i].c.x);;
		fprintf(fp2, "%f ",list.item[i].a.y);	
		fprintf(fp2, "%f ",list.item[i].b.y);	
		fprintf(fp2, "%f\n",list.item[i].c.y);
		count = 0;
		for(j=0; j<np; j++){
			if((p[j].x==list.item[i].a.x)&&(p[j].y==list.item[i].a.y)){
				fprintf(fp3,"%d ",j);
				count++;
			}
		}
		for(j=0; j<np; j++){
			if((p[j].x==list.item[i].b.x)&&(p[j].y==list.item[i].b.y)){
				fprintf(fp3,"%d ",j);
				count++;
			}
		}
		for(j=0; j<np; j++){
			if((p[j].x==list.item[i].c.x)&&(p[j].y==list.item[i].c.y)){
				fprintf(fp3,"%d\n",j);
				count++;
			}
		}
	}
	fclose(fp1);
	fclose(fp2);
	free(p);
	tlistfree(&list);
}

