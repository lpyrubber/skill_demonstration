#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

#define gamma 1.4
#define CFL 0.5
#define Mach 2.5
#define P_0 10000.0
#define R 287.0
#define T_0 300.0
#define CV (R/(gamma-1))
#define CP (CV+R)

#define STEP 1000
#define TOR 1e-9
#define N_it 10000

class cell{
public:
    int slist[4];
    int nbr[4];
    Eigen::Matrix<double, 4, 1> U, dU, dU_old, Res;
    Eigen::Matrix<double, 4, 4> M, Mg;
    double x, y, vol;
    double rho, P, u, v, a, T; 
};

class surface{
public:
    double nx, ny, x, y, l;
    int nbr[2];
    Eigen::Matrix<double, 4, 4> Am, Ap;
    Eigen::Matrix<double, 4, 1> F;
    int label;
};

void Load_File();
void Allocate_Memory();
void Initial();
void Calculate_Flux();
void Interface_Calculation(int surf);
void Update_U();
void Save_Result();
void Free_Memory();


int Nx, Ny, Nc, Np, Ns;
double **xpoints;
double *res;
double dt;
double temp=0;
double Tt=0;
Eigen::Matrix<double, 4, 4> E;
Eigen::Matrix<double, 4, 1> x, b;
Eigen::Matrix<double, 4, 1> row;
Eigen::Matrix<double, 1, 4> col;
Eigen::Matrix<double, 4, 1> Ug;
cell* Cell;
surface* Surface;

int main() {
    int i;
    FILE* in;
    in = fopen("residual.txt", "w");
    Load_File();
    dt=(double)(CFL*1.0/(Mach*Nx*400));
    Tt=dt;
    Initial();
    
    printf("dt=%f\n",dt);
    for(i=0; i<STEP; i++){
        Calculate_Flux();
        Update_U();
        fprintf(in, "%le %le\n", Tt, temp);
        
    }
    fprintf(in, "%le %le\n", Tt, temp);
    fclose(in);
    Save_Result();
    Free_Memory();
    return 0;
}


void Load_File(){
    FILE* in;
    FILE* in2;
    FILE* in3;
    int i,j;
    int temp1, temp2;
    in = fopen("grid_info.txt", "r");
    fscanf(in, "%d %d %d %d %d\n",&Np, &Nx, &Ny, &Nc, &Ns);
    fclose(in);
    Allocate_Memory();
    printf("%d %d %d %d %d\n",Np, Nx, Ny, Nc, Ns);

    in = fopen("cell_info.txt","r");
    in2 = fopen("clist.txt","r");
    for(i=0; i<Nc; i++){
        fscanf(in, "%le %le %le\n", &Cell[i].vol, &Cell[i].x, &Cell[i].y);
        fscanf(in2, "%d %d %d %d %d %d %d %d\n", Cell[i].slist, Cell[i].slist+1, Cell[i].slist+2, Cell[i].slist+3, Cell[i].nbr, Cell[i].nbr+1, Cell[i].nbr+2, Cell[i].nbr+3);

    }
    fclose(in);
    fclose(in2);


    in = fopen("surface_info.txt","r");
    in2 = fopen("slist.txt","r");    
    for(i=0; i<Ns; i++){
        fscanf(in, "%le %le %le %le %le\n", &Surface[i].l, &Surface[i].nx, &Surface[i].ny, &Surface[i].x, &Surface[i].y);
        fscanf(in2, "%d %d %d %d %d\n", &temp1, &temp2, Surface[i].nbr, Surface[i].nbr+1, &Surface[i].label);
    }

    fclose(in);
    fclose(in2);


}

void Allocate_Memory(){
    int i;
    Cell = new cell[Nc];
    Surface = new surface[Ns];
}

void Initial(){
    int i, j;
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){

            Cell[i+j*Nx].T=T_0;
            Cell[i+j*Nx].a=sqrt(gamma*R*T_0);
            Cell[i+j*Nx].u=Mach*Cell[i+j*Nx].a*j/Ny;
            Cell[i+j*Nx].v=0;
            Cell[i+j*Nx].P=P_0;
            Cell[i+j*Nx].rho=P_0/(R*T_0);
 
            Cell[i+j*Nx].U(0,0)=Cell[i+j*Nx].rho;
            Cell[i+j*Nx].U(1,0)=Cell[i+j*Nx].rho*Cell[i+j*Nx].u;
            Cell[i+j*Nx].U(2,0)=Cell[i+j*Nx].rho*Cell[i+j*Nx].v;
            Cell[i+j*Nx].U(3,0)=Cell[i+j*Nx].rho*(CV*Cell[i+j*Nx].T+0.5*(Cell[i+j*Nx].u*Cell[i+j*Nx].u+Cell[i+j*Nx].v*Cell[i+j*Nx].v));
            //Cell[i+j*Nx].Nbr[0]=0;
            Cell[i+j*Nx].Res.setZero();
            Cell[i+j*Nx].dU_old.setZero();
            Cell[i+j*Nx].dU=Cell[i+j*Nx].dU_old;
            Cell[i+j*Nx].M=Cell[i+j*Nx].vol/dt*Cell[i+j*Nx].M.setIdentity();
//            Cell[i+j*Nx].Mg=Cell[i+j*Nx].M;
        }
    }
}

void Interface_Calculation(int surf){

    double u1, u2, u3, u4, u, v, p, a, nx, ny, ud;
    double lt, lc, Prho, PE, h0, la, lap, lam;
    double u1g, u2g, u3g, u4g, rhog, ug, vg, pg;
    int i,j, Lid, Rid;
    Lid=Surface[surf].nbr[0];
    Rid=Surface[surf].nbr[1];
    if(Surface[surf].label==0){
        u1 = 0.5*(Cell[Lid].U(0,0)+Cell[Rid].U(0,0));
        u2 = 0.5*(Cell[Lid].U(1,0)+Cell[Rid].U(1,0));
        u3 = 0.5*(Cell[Lid].U(2,0)+Cell[Rid].U(2,0));
        u4 = 0.5*(Cell[Lid].U(3,0)+Cell[Rid].U(3,0));
        nx=Surface[surf].nx;
        ny=Surface[surf].ny;

    }else{
        if((Surface[surf].label==1)||(Surface[surf].label==3)){
            ug=Cell[Lid].u-2*(Cell[Lid].u*Surface[surf].nx+Cell[Lid].v*Surface[surf].ny)*Surface[surf].nx;
            vg=Cell[Lid].v-2*(Cell[Lid].u*Surface[surf].nx+Cell[Lid].v*Surface[surf].ny)*Surface[surf].ny;
            rhog=Cell[Lid].rho;
            Ug(0,0)=rhog;
            Ug(1,0)=rhog*ug;
            Ug(2,0)=rhog*vg;
            Ug(3,0)=Cell[Lid].U(3);
            
        }else if(Surface[surf].label==2){
            Ug = Cell[Lid].U;
        }else{
            rhog=P_0/(R*T_0);
            ug=Mach*sqrt(gamma*R*T_0);
            Ug(0,0)=rhog;
            Ug(1,0)=rhog*ug;
            Ug(2,0)=0;
            Ug(3,0)=rhog*(CV*T_0+0.5*(ug*ug));
        }
        u1 = 0.5*(Cell[Lid].U(0,0)+Ug(0,0));
        u2 = 0.5*(Cell[Lid].U(1,0)+Ug(1,0));
        u3 = 0.5*(Cell[Lid].U(2,0)+Ug(2,0));
        u4 = 0.5*(Cell[Lid].U(3,0)+Ug(3,0));
        nx=Surface[surf].nx;
        ny=Surface[surf].ny;
    }

    u = u2/u1;
    v = u3/u1;
    p = (u4-0.5*u1*(u*u+v*v))*(gamma-1);
    a = sqrt(gamma*p/u1);

    ud=u*nx+v*ny;   
    Prho=0.5*(gamma-1)*(u*u+v*v);
    PE=(gamma-1);
    h0=(u4+p)/u1;
    //Ap
    la=(ud>0)?ud:0;
    lap=((ud+a)>0)?ud+a:0;
    lam=((ud-a)>0)?ud-a:0;
    lt=0.5*(lap-lam);
    lc=0.5*(lap+lam)-la;

    row(0,0)=lc/(a*a);
    row(1,0)=(u*lc+a*nx*lt)/(a*a);
    row(2,0)=(v*lc+a*ny*lt)/(a*a);
    row(3,0)=(h0*lc+a*ud*lt)/(a*a);
    col(0,0)=Prho;
    col(0,1)=-u*PE;
    col(0,2)=-v*PE;
    col(0,3)=PE;
    Surface[surf].Ap=row*col;

    row(0,0)=lt/a;
    row(1,0)=(u*lt/a)+nx*lc;
    row(2,0)=(v*lt/a)+ny*lc;
    row(3,0)=(h0*lt/a)+ud*lc;
    col(0,0)=-ud;
    col(0,1)=nx;
    col(0,2)=ny;
    col(0,3)=0;
    Surface[surf].Ap+=row*col;
    for(j=0; j<4; j++){
        Surface[surf].Ap(j,j)+=la;
    }

    //Am
    la=(ud<0)?ud:0;
    lap=((ud+a)<0)?ud+a:0;
    lam=((ud-a)<0)?ud-a:0;
    lt=0.5*(lap-lam);
    lc=0.5*(lap+lam)-la;
    row(0,0)=lc/(a*a);
    row(1,0)=(u*lc+a*nx*lt)/(a*a);
    row(2,0)=(v*lc+a*ny*lt)/(a*a);
    row(3,0)=(h0*lc+a*ud*lt)/(a*a);
    col(0,0)=Prho;
    col(0,1)=-u*PE;
    col(0,2)=-v*PE;
    col(0,3)=PE;
    Surface[surf].Am=row*col;

    row(0,0)=lt/a;
    row(1,0)=(u*lt/a)+nx*lc;
    row(2,0)=(v*lt/a)+ny*lc;
    row(3,0)=(h0*lt/a)+ud*lc;
    col(0,0)=-ud;
    col(0,1)=nx;
    col(0,2)=ny;
    col(0,3)=0;
    Surface[surf].Am+=row*col;
    for(j=0; j<4; j++){
        Surface[surf].Am(j,j)+=la;
    }
    Surface[surf].Ap=Surface[surf].l*Surface[surf].Ap;
    Surface[surf].Am=Surface[surf].l*Surface[surf].Am;
    x=Surface[surf].Ap*Cell[Lid].U;
    if(Surface[surf].label==0){
        x+=Surface[surf].Am*Cell[Rid].U;
        Surface[surf].F=x;
        Cell[Lid].M=Cell[Lid].M+Surface[surf].Ap;
 //       Cell[Lid].Mg=Cell[Lid].Mg+Surface[surf].Ap;
        Cell[Lid].Res=Cell[Lid].Res-x;
//        Cell[Lid].N[Cell[Lid].Nbr[0]]=Surface[surf].Am;
//        Cell[Lid].Nbr[0]++;
//        Cell[Lid].Nbr[Cell[Lid].Nbr[0]]=Rid;

        Cell[Rid].M=Cell[Rid].M-Surface[surf].Am;
 //       Cell[Rid].M=Cell[Rid].Mg-Surface[surf].Am;
        Cell[Rid].Res=Cell[Rid].Res+x;
//        Cell[Rid].N[Cell[Rid].Nbr[0]]=-1*Surface[surf].Ap;
//        Cell[Rid].Nbr[0]++;
//        Cell[Rid].Nbr[Cell[Rid].Nbr[0]]=Lid;


    }else{
        x+=Surface[surf].Am*Ug;      
        E.setIdentity();
        if(Surface[surf].label==1||Surface[surf].label==3){
            E(1,1)-=2*nx*nx;
            E(1,2)-=2*nx*ny;
            E(2,1)-=2*ny*nx;
            E(2,2)-=2*ny*ny;
        }else if(Surface[surf].label==4){
            E.setZero();
        }
        //supersonic outflow, E=I
        Cell[Lid].M=Cell[Lid].M+Surface[surf].Ap+Surface[surf].Am*E;
  //      Cell[Lid].Mg=Cell[Lid].Mg+Surface[surf].Ap;
        Surface[surf].F=x;
        Cell[Lid].Res=Cell[Lid].Res-x;

    }

}

void Calculate_Flux(){
    int j;
    for(j=0; j<Ns; j++){
        Interface_Calculation(j);
    }
//    printf("Finish Flux calculation\n");
}

void Update_U(){
    int i=0, j, k, surf;


     for(j=0; j<Nc; j++){
        for(k=0; k<4; k++){
            surf=Cell[j].slist[k];
            if(j==Surface[surf].nbr[0]){
                Cell[j].U-=dt/Cell[j].vol*Surface[surf].F;
            }else{
                Cell[j].U+=dt/Cell[j].vol*Surface[surf].F;
            }
//          b-=Cell[j].N[k]*Cell[Cell[j].Nbr[k+1]].dU_old;
        }
    }
    i++;
 /*   
    for(j=0; j<Nc; j++){
        Cell[j].dU_old=Cell[j].Mg.inverse()*Cell[j].Res;
    }
 */

   /*
    while((residual>TOR)&&(i<N_it)){
        residual=0;
        for(j=0; j<Nc; j++){
            b = Cell[j].Res;
            for(k=0; k<4; k++){
                surf=Cell[j].slist[k];
                if(Surface[surf].label==0){
                    if(j==Surface[surf].nbr[0]){
                        b-=Surface[surf].Am*Cell[Surface[surf].nbr[1]].dU_old;
                    }else{
                        b+=Surface[surf].Ap*Cell[Surface[surf].nbr[0]].dU_old;
                    }
                }
//                b-=Cell[j].N[k]*Cell[Cell[j].Nbr[k+1]].dU_old;
            }

            //std::cout<<b<<std::endl;
            Cell[j].dU=Cell[j].M.inverse()*b;
            residual+=fabs(Cell[j].dU(0,0)-Cell[j].dU_old(0,0))\
                     +fabs(Cell[j].dU(1,0)-Cell[j].dU_old(1,0))\
                     +fabs(Cell[j].dU(2,0)-Cell[j].dU_old(2,0))\
                     +fabs(Cell[j].dU(3,0)-Cell[j].dU_old(3,0));
            Cell[j].dU_old=Cell[j].dU;

        }
        i++;
    }
    */
    //printf("residual=%f\n",residual);
    double RES=0;
    temp=0;
    double max=0;
    for(i=0; i<Nc; i++){
        temp+=Cell[i].Res(0)*Cell[i].Res(0)/Cell[i].vol*Cell[i].vol;
        Cell[i].rho=Cell[i].U(0,0);
        Cell[i].u=Cell[i].U(1,0)/Cell[i].rho;
        Cell[i].v=Cell[i].U(2,0)/Cell[i].rho;
        Cell[i].T=(Cell[i].U(3,0)/Cell[i].rho-0.5*(Cell[i].u*Cell[i].u+Cell[i].v*Cell[i].v))/CV;
        Cell[i].P=Cell[i].rho*R*Cell[i].T;
        Cell[i].a=sqrt(gamma*R*Cell[i].T);
        Cell[i].dU.setZero();
        Cell[i].dU_old=Cell[i].dU;
        Cell[i].Res.setZero();
        Cell[i].M=Cell[i].vol/dt*Cell[i].M.setIdentity();
        Cell[i].Mg=Cell[i].M;
        max=((fabs(Cell[i].u)+Cell[i].a)>max)?(fabs(Cell[i].u)+Cell[i].a):max;
    }   
    dt=(double)(CFL)/(max*double(Nx));
    Tt+=dt;
    printf("max=%f, res=%f dt=%f\n", max, sqrt(temp), dt);

//    printf("Finish Update U\n");
}

void Save_Result(){
    int i,j;
    double m;
    FILE* in;
    in = fopen("result.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            fprintf(in, "%le %le %le %le\n",Cell[i+j*Nx].rho, Cell[i+j*Nx].T, Cell[i+j*Nx].P, Cell[i+j*Nx].U(3,0));
        }
    }
    fclose(in);

    in = fopen("mesh.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            fprintf(in, "%le %le\n",Cell[i+j*Nx].x, Cell[i+j*Nx].y);
        }
    }
    fclose(in);

    in = fopen("u.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            m=sqrt(Cell[i+j*Nx].u * Cell[i+j*Nx].u + Cell[i+j*Nx].v * Cell[i+j*Nx].v)/Cell[i+j*Nx].a;
            fprintf(in, "%le %le %le\n",Cell[i+j*Nx].u, Cell[i+j*Nx].v,m);
        }
    }
    fclose(in);


}

void Free_Memory(){
    int i;
    delete[] Cell;
    delete[] Surface;
}

/*
    x=Eigen::Matrix<double,4,1>::Zero();
    //poesitive nx
    int ind, nbr,surf;

    ind=101;
    nbr=Cell[ind].nbr[1];
    surf=Cell[ind].slist[1];


    //negative nx
    ind=102;
    nbr=Cell[ind].nbr[0];
    surf=Cell[ind].slist[0];


    
    Cell[ind].dU_old(0,0)=0.266;
    Cell[ind].dU_old(1,0)=0.266*0.927;
    Cell[ind].dU_old(2,0)=0;

    Cell[ind].dU_old(3,0)=0.303/(gamma-1)+0.5*0.266*0.927*0.927;

    printf("L: %f %f %f %f\n",Cell[ind].dU_old(0,0), Cell[ind].dU_old(1,0), Cell[ind].dU_old(2,0), Cell[ind].dU_old(3,0));


    Cell[nbr].dU_old(0,0)=0.125;
    Cell[nbr].dU_old(1,0)=0;
    Cell[nbr].dU_old(2,0)=0;
    Cell[nbr].dU_old(3,0)=0.1/(gamma-1);

    
    printf("L: %f %f %f %f\n",Cell[nbr].dU_old(0,0), Cell[nbr].dU_old(1,0), Cell[nbr].dU_old(2,0), Cell[nbr].dU_old(3,0));
    
    printf("dist:%f %f\n",Cell[nbr].x-Cell[ind].x,Cell[nbr].y-Cell[ind].y);
    printf("n:%f %f\n",Surface[surf].nx, Surface[surf].ny);
    printf("surf nbr: %d %d\n",Surface[surf].nbr[0],Surface[surf].nbr[1]);

    Interface_Calculation(ind, nbr, surf);
*/
