#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NU    0.3    //Poisson's ratio
#define E     21e9   //Modulus of Elasticity
#define THICK 0.025  //Thickness of plate

void Calc_Local_K(float *k, float x1, float y1, float x2, float y2, float x3, float y3);
void Calc_Global_K(float *K, float *k, int ele1, int ele2, int ele3);
void Partition_K_F();
void Calc_Stress*();

int main(){
	float *x, *y, *Fx, *Fy, *Fix_x, *Fix_y, *tri, *A, *B, *D, *K, *k, *Kr, *Fr;
	int np, ne;
	//get np, ne;
	char buffer[256];
	FILE *fp;

	
