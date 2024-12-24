#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ANGLE 0 //in degree
#define LX 1.0   //inlet
#define LY 0.25  //length
#define SX (1.0/3.0)  //slope start point
#define RX (1.0/3.0)  //x-direction length for the slope
#define NX 50 
#define NY 25
#define DX (double)(LX/NX)
#define DY (double)(LY/NY)
#define DIM 2

void Allocate_Memory();
void Generate_Grid();
void Calculate_Surface();
void Calculate_Cell();
void Save_Result();
void Free_Memory();

double *xpoint, *normal_vector, *length, *xcell, *volume, *skew;  
float t_max, t_min;
int **clist, **slist, **plist, *slabel;
int Np, Nc, Ns;


int main(){
    Np=(NX+1)*(NY+1);
    Nc=NX*NY;
    Ns=NX*(NY+1)+(NX+1)*NY;
    Allocate_Memory();
    Generate_Grid();
    Calculate_Surface();
    Calculate_Cell();
    Save_Result();
    Free_Memory();
}


void Allocate_Memory(){
    int i,j;
    //double
    xpoint=(double*)malloc(Np*DIM*sizeof(double));
    volume=(double*)malloc(Nc*sizeof(double));
    xcell=(double*)malloc(Nc*DIM*sizeof(double));
    normal_vector=(double*)malloc(2*Ns*DIM*sizeof(double));
    length=(double*)malloc(Ns*sizeof(double));
    skew=(double*)malloc(Np*sizeof(double));
    
    //int
    clist=(int**)malloc(Nc*sizeof(int*)); //4 faces
    plist=(int**)malloc(Nc*sizeof(int*)); //4 faces
    slist=(int**)malloc(Ns*sizeof(int*)); //2 point
    slabel=(int*)malloc(Ns*sizeof(int));
    for(i=0; i<Nc; i++){
        clist[i]=(int*)malloc(8*sizeof(int));
        plist[i]=(int*)malloc(4*sizeof(int));
    }
    for(i=0; i<Ns; i++){
        slist[i]=(int*)malloc(4*sizeof(int));
    }
    
}

void Generate_Grid(){
    int i, j;
    double dyr, dys, ur, dr, us, ds;
    us=LY-(RX*tan(ANGLE*M_PI/180));
    ds=LY-us;
    dys=us/NY;
    for(i=0; i<(NX+1); i++){
        ur=LY-((i*DX-SX)*tan(ANGLE*M_PI/180));
        dr=LY-ur;
        dyr=ur/NY;
        for(j=0; j<(NY+1); j++){
            xpoint[i+j*(NX+1)]=i*DX;
            if(i*DX<SX){
                xpoint[i+j*(NX+1)+Np]=j*DY;
            }else if(i*DX<(SX+RX)){
                xpoint[i+j*(NX+1)+Np]=j*dyr+dr;
            }else{
                xpoint[i+j*(NX+1)+Np]=j*dys+ds;
            }
        }
    }

    //update plist
    for(j=0; j<NY; j++){
        for(i=0; i<NX; i++){
            //counter clockwise [0,0]->[1,0]->[1,1]->[0,1]
            plist[i+j*NX][0]=i+j*(NX+1);
            plist[i+j*NX][1]=i+1+j*(NX+1);
            plist[i+j*NX][2]=i+1+(j+1)*(NX+1);
            plist[i+j*NX][3]=i+(j+1)*(NX+1);
        }
    }
}


void Calculate_Surface(){
    int index=0;
    int i,j;
    //create surface list
    //y-dir
    for(j=0; j<NY; j++){
        //x-dir
        for(i=0; i<NX; i++){
            //right node
            slist[index][0]=i+j*(NX+1);
            slist[index][1]=i+1+j*(NX+1);
            slabel[index]=(j==0)?1:0;
            index++;
            
            //up node
            slist[index][0]=i+j*(NX+1);
            slist[index][1]=i+(j+1)*(NX+1);
            slabel[index]=(i==0)?4:0;
            index++;
        }
        //up edge
        slist[index][0]=NX+j*(NX+1);
        slist[index][1]=NX+(j+1)*(NX+1);
        slabel[index]=2;
        index++;
    }
    //floor edge
    for(i=0; i<NX; i++){
        slist[index][0]=i+NY*(NX+1);
        slist[index][1]=i+1+NY*(NX+1);
        slabel[index]=3;
        index++;
    }
    //calculate length and normal vector
    for(i=0; i<Ns; i++){
        length[i]=sqrt((xpoint[slist[i][0]]-xpoint[slist[i][1]])*(xpoint[slist[i][0]]-xpoint[slist[i][1]])+(xpoint[slist[i][0]+Np]-xpoint[slist[i][1]+Np])*(xpoint[slist[i][0]+Np]-xpoint[slist[i][1]+Np]));
    }
    //2d normal: (x,y)->(y,-x)
    for(i=0; i<Ns; i++){
        normal_vector[i]=(xpoint[slist[i][0]+Np]-xpoint[slist[i][1]+Np])/length[i];
        normal_vector[i+Ns]=-(xpoint[slist[i][0]]-xpoint[slist[i][1]])/length[i];
        normal_vector[i+2*Ns]=(xpoint[slist[i][0]]+xpoint[slist[i][1]])*0.5;
        normal_vector[i+3*Ns]=(xpoint[slist[i][0]+Np]+xpoint[slist[i][1]+Np])*0.5;
    }
}

void Calculate_Cell(){
    int i,j,k, ix, iy, j1, j2, j3, sindex, max_id=-1, min_id=-1;
    double cx, cy, temp, temp1, temp2, nbrx, nbry, delta, x1, x2, y1, y2;
    //update surface to cell
    for(i=0; i<Nc; i++){
        for(j=0; j<4; j++){
            j2=(j+1)%4;
            if(plist[i][j]>plist[i][j2]){
                ix=plist[i][j2];
                iy=plist[i][j];
            }else{
                ix=plist[i][j];
                iy=plist[i][j2];
            }
            k=0;
            while(slist[k][0]!=ix){
                k++;
            }
            while(slist[k][1]!=iy){
                k++;
            }
            clist[i][j]=k;
        }
    }
    //calculate x cell center
    for(i=0; i<Nc; i++){
        cx=0;
        cy=0;
        for(j=0; j<4; j++){
            cx+=xpoint[plist[i][j]];
            cy+=xpoint[plist[i][j]+Np];
        }
        xcell[i]=0.25*cx;
        xcell[i+Nc]=0.25*cy;
    }
    for(i=0; i<Nc; i++){
        temp=0;
        for(j=0; j<2; j++){
            
            j1=2*j;
            j2=(j1+1)%4;
            j3=(j1-1+4)%4;
            //0.5(x1y2-y1x2)
            temp+=0.5*fabs((xpoint[plist[i][j2]]-xpoint[plist[i][j1]])*(xpoint[plist[i][j3]+Np]-xpoint[plist[i][j1]+Np])-(xpoint[plist[i][j2]+Np]-xpoint[plist[i][j1]+Np])*(xpoint[plist[i][j3]]-xpoint[plist[i][j1]]));
        }
        volume[i]=temp;
    }
    t_max=0;
    t_min=1e6;
    //question 1 required value
    for(j=0; j<NY; j++){
        for(i=0; i<NX; i++){
            temp1=0;
            temp2=0;
            sindex=clist[i+j*NX][0];
            if(j==0){
                //mirror the cell center with the boundary surface
                x1=xpoint[slist[sindex][0]];
                y1=xpoint[slist[sindex][0]+Np];
                x2=xpoint[slist[sindex][1]];
                y2=xpoint[slist[sindex][1]+Np];
                nbrx=2*(x1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(x2-x1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX];
                nbry=2*(y1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(y2-y1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX+Nc];
                clist[i+j*NX][4]=-slabel[sindex];
            }else{
                nbrx=xcell[i+(j-1)*NX];
                nbry=xcell[i+(j-1)*NX+Nc];
                clist[i+j*NX][4]=i+(j-1)*NX;
            }
            //(xnbr-xcell)*normalvector
            delta=(nbrx-xcell[i+j*NX])*normal_vector[sindex]+(nbry-xcell[i+j*NX+Nc])*normal_vector[sindex+Ns];
            delta=delta/fabs(delta);

            
            
            temp1+=delta*normal_vector[sindex]*length[sindex];
            temp2+=delta*normal_vector[sindex+Ns]*length[sindex];

            if(slabel[sindex]){
                    normal_vector[sindex]=delta*normal_vector[sindex];
                    normal_vector[sindex+Ns]=delta*normal_vector[sindex+Ns];
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=-1;
                
            }else{
                if(delta>0){
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=clist[i+j*NX][4];
                }
            }

            sindex=clist[i+NX*j][1];
            if(i==NX-1){
                //mirror the cell center with the boundary surface
                x1=xpoint[slist[sindex][0]];
                y1=xpoint[slist[sindex][0]+Np];
                x2=xpoint[slist[sindex][1]];
                y2=xpoint[slist[sindex][1]+Np];
                nbrx=2*(x1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(x2-x1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX];
                nbry=2*(y1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(y2-y1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX+Nc];
                clist[i+j*NX][5]=-slabel[sindex];
            }else{
                nbrx=xcell[(i+1)+j*NX];
                nbry=xcell[(i+1)+j*NX+Nc];
                clist[i+j*NX][5]=(i+1)+j*NX;
            }
            //(xnbr-xcell)*normalvector
            delta=(nbrx-xcell[i+j*NX])*normal_vector[sindex]+(nbry-xcell[i+j*NX+Nc])*normal_vector[sindex+Ns];
            delta=delta/fabs(delta);
            

            temp1+=delta*normal_vector[sindex]*length[sindex];
            temp2+=delta*normal_vector[sindex+Ns]*length[sindex];

            if(slabel[sindex]){
                    normal_vector[sindex]=delta*normal_vector[sindex];
                    normal_vector[sindex+Ns]=delta*normal_vector[sindex+Ns];
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=-1;
                
            }else{
                if(delta>0){
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=clist[i+j*NX][5];
                }
            }

            sindex=clist[i+NX*j][2];
            if(j==NY-1){
                //mirror the cell center with the boundary surface
                x1=xpoint[slist[sindex][0]];
                y1=xpoint[slist[sindex][0]+Np];
                x2=xpoint[slist[sindex][1]];
                y2=xpoint[slist[sindex][1]+Np];
                nbrx=2*(x1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(x2-x1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX];
                nbry=2*(y1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(y2-y1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX+Nc];
                clist[i+j*NX][6]=-slabel[sindex];
            }else{
                nbrx=xcell[i+(j+1)*NX];
                nbry=xcell[i+(j+1)*NX+Nc];
                clist[i+j*NX][6]=i+(j+1)*NX;
            }
            //(xnbr-xcell)*normalvector
            delta=(nbrx-xcell[i+j*NX])*normal_vector[sindex]+(nbry-xcell[i+j*NX+Nc])*normal_vector[sindex+Ns];
            delta=delta/fabs(delta);
            
            temp1+=delta*normal_vector[sindex]*length[sindex];
            temp2+=delta*normal_vector[sindex+Ns]*length[sindex];

            if(slabel[sindex]){
                    normal_vector[sindex]=delta*normal_vector[sindex];
                    normal_vector[sindex+Ns]=delta*normal_vector[sindex+Ns];
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=-1;
                
            }else{
                if(delta>0){
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=clist[i+j*NX][6];
                }
            }

            sindex=clist[i+NX*j][3];
            if(i==0){
                //mirror the cell center with the boundary surface
                x1=xpoint[slist[sindex][0]];
                y1=xpoint[slist[sindex][0]+Np];
                x2=xpoint[slist[sindex][1]];
                y2=xpoint[slist[sindex][1]+Np];
                nbrx=2*(x1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(x2-x1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX];
                nbry=2*(y1+((xcell[i+j*NX]-x1)*(x2-x1)+(xcell[i+j*NX+Nc]-y1)*(y2-y1))*(y2-y1)/(pow(x2-x1,2)+pow(y2-y1,2)))-xcell[i+j*NX+Nc];
                clist[i+j*NX][7]=-slabel[sindex];
            }else{
                nbrx=xcell[(i-1)+j*NX];
                nbry=xcell[(i-1)+j*NX+Nc];
                clist[i+j*NX][7]=(i-1)+j*NX;
            }
            //(xnbr-xcell)*normalvector
            delta=(nbrx-xcell[i+j*NX])*normal_vector[sindex]+(nbry-xcell[i+j*NX+Nc])*normal_vector[sindex+Ns];
            delta=delta/fabs(delta);
            
            temp1+=delta*normal_vector[sindex]*length[sindex];
            temp2+=delta*normal_vector[sindex+Ns]*length[sindex];

            if(slabel[sindex]){
                    normal_vector[sindex]=delta*normal_vector[sindex];
                    normal_vector[sindex+Ns]=delta*normal_vector[sindex+Ns];
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=-1;
                
            }else{
                if(delta>0){
                    slist[sindex][2]=i+j*NX;
                    slist[sindex][3]=clist[i+j*NX][7];
                }
            }

          
            temp=sqrt(temp1*temp1+temp2*temp2);
            skew[i+j*NX]=temp;
            if(temp>t_max){
                t_max=temp;
                max_id=i+j*NX;
            }
            if(temp<t_min){
                t_min=temp;
                min_id=i+j*NX;
            }
        }
    }
    printf("max=%e at (%d,%d), min=%e at (%d,%d)\n", t_max, max_id%NX, max_id/NX, t_min, min_id%NX, min_id/NX);
}

void Save_Result(){

    FILE* out;
    int i,j;
    double sum;
    sum=0;
    out = fopen("ramp_points.txt","w");
    for(i=0; i<Np; i++){
        fprintf(out, "%e %e\n", xpoint[i], xpoint[i+Np]);
    }
    fclose(out);

    out = fopen("plist.txt","w");
    for(i=0; i<Nc; i++){
        fprintf(out, "%d %d %d %d\n", plist[i][0], plist[i][1], plist[i][2], plist[i][3]);
    }
    fclose(out);

    out = fopen("clist.txt","w");
    for(i=0; i<Nc; i++){
        fprintf(out, "%d %d %d %d %d %d %d %d\n", clist[i][0], clist[i][1], clist[i][2], clist[i][3], clist[i][4], clist[i][5], clist[i][6], clist[i][7]);
    }
    fclose(out);

    out = fopen("cell_info.txt","w");
    for(i=0; i<Nc; i++){
        sum+=volume[i];
        fprintf(out, "%e %e %e\n", volume[i], xcell[i], xcell[i+Nc]);
    }
    printf("total volume=%e\n",sum);
    fclose(out);

    out = fopen("slist.txt","w");
    for(i=0; i<Ns; i++){
        fprintf(out, "%d %d %d %d %d\n", slist[i][0], slist[i][1], slist[i][2], slist[i][3],slabel[i]);
    }
    fclose(out);

    out = fopen("surface_info.txt","w");
    for(i=0; i<Ns; i++){
        fprintf(out, "%e %e %e %e %e\n", length[i], normal_vector[i], normal_vector[i+Ns], normal_vector[i+2*Ns], normal_vector[i+3*Ns]);
    }
    fclose(out);

    out = fopen("grid_info.txt","w");
    fprintf(out, "%d %d %d %d %d\n", Np, NX, NY, Nc, Ns);
}

void Free_Memory(){
    int i;
    free(xpoint);
    free(volume);
    free(xcell);
    free(normal_vector);
    free(length);
    for(i=0; i<Nc; i++){
        free(clist[i]);
    }
    free(clist);
    for(i=0; i<Ns; i++){
        free(slist[i]);
    }
    free(slist);
    free(slabel);

}
