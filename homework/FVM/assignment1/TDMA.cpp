#include <stdlib.h>
#include <stdio.h>

#define N 41
#define L 1.0
#define NIT 1000

float *x,*d,*cp,*dp;

int main(){
	int i,j,k;
	float c,b,a;
	x = (float*)malloc(N*sizeof(float));
	d = (float*)malloc(N*sizeof(float));
	cp = (float*)malloc(N*sizeof(float));
	dp = (float*)malloc(N*sizeof(float));
	x[0]=0;
	x[N-1]=1;
	c=1.0;
	b=-2.0;
	a=1.0;
	for(i=1;i<N-2;i++){
		d[i]=0;		
	}
	d[N-2]=-1;
	cp[1]=c/b;
	dp[1]=d[1]/b;
	for(i=2; i<N-1; i++){
		cp[i]=c/(b-a*cp[i-1]);
		dp[i]=(d[i]-a*dp[i-1])/(b-a*cp[i-1]);	
	}
	x[N-2]=dp[N-2];
	for(i=N-3;i>0;i--){
		x[i]=dp[i]-cp[i]*x[i+1];
	}
	for(i=0; i<N; i++){
		printf("%f ",x[i]);
	}
	free(x);
	free(d);
	free(cp);
	free(dp);
	return 0;
}
