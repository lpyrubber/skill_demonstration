#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <time.h>
#include <omp.h>
/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif


#define gamma 1.4
#define CFL 4
#define Mach 2.5
#define P_0 10000.0
#define R 287.0
#define T_0 300.0
#define CV (R/(gamma-1))
#define CP (CV+R)


#define STEP 80

class cell{
public:
    int slist[4];
    int nbr[4];
    Eigen::Matrix<double, 4, 1> U, Res;
    Eigen::Matrix<double, 4, 4> M;
    double x, y, vol;
    double rho, P, u, v, a, T; 
};

class surface{
public:
    double nx, ny, x, y, l;
    int nbr[2];
    int label;
};

void Load_File();
void Allocate_Memory();
void Initial();
void Calculate_Flux();
void Interface_Calculation(int Lid, int Rid, int surf);
void Gather_Matrix();
void Update_U();
void Save_Result();
void Free_Memory();


int Nx, Ny, Nc, Np, Ns;
double **xpoints;
double *res, dt;
Eigen::SparseMatrix<double> A(10000,10000);
Eigen::VectorXd xi(10000), bi(10000);
Eigen::Matrix<double, 4, 4> E, Ap, Am;
Eigen::Matrix<double, 4, 1> x, b;
Eigen::Matrix<double, 4, 1> row;
Eigen::Matrix<double, 1, 4> col;
Eigen::Matrix<double, 4, 1> Ug;
cell* Cell;
surface* Surface;

int main() {


    int i;
    Load_File();
    dt=(double)(CFL*1.0/(Mach*Nx*400));
    Initial();
  
    for(i=0; i<STEP; i++){
        Calculate_Flux();
        Update_U();
        printf("finish %d step\n",i);
    }

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
    A.resize(Nc*4,Nc*4);
    A.reserve(Eigen::VectorXi::Constant(Nc*4,5*4));
    xi.resize(Nc*4);
    bi.resize(Nc*4);
    Cell = new cell[Nc];
    Surface = new surface[Ns];
}

void Initial(){
    int i, j;
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){

            Cell[i+j*Nx].T=T_0;
            Cell[i+j*Nx].a=sqrt(gamma*R*T_0);
            Cell[i+j*Nx].u=Mach*Cell[i+j*Nx].a;
            Cell[i+j*Nx].v=0;
            Cell[i+j*Nx].P=P_0;
            Cell[i+j*Nx].rho=P_0/(R*T_0);
 
            Cell[i+j*Nx].U(0,0)=Cell[i+j*Nx].rho;
            Cell[i+j*Nx].U(1,0)=Cell[i+j*Nx].rho*Cell[i+j*Nx].u;
            Cell[i+j*Nx].U(2,0)=Cell[i+j*Nx].rho*Cell[i+j*Nx].v;
            Cell[i+j*Nx].U(3,0)=Cell[i+j*Nx].rho*(CV*Cell[i+j*Nx].T+0.5*(Cell[i+j*Nx].u*Cell[i+j*Nx].u+Cell[i+j*Nx].v*Cell[i+j*Nx].v));
            //Cell[i+j*Nx].Nbr[0]=0;
            Cell[i+j*Nx].Res.setZero();
            Cell[i+j*Nx].M=Cell[i+j*Nx].vol/dt*Cell[i+j*Nx].M.setIdentity();
        }
    }
}

void Interface_Calculation(int surf){

    double u1, u2, u3, u4, u, v, p, a, nx, ny, ud;
    double lt, lc, Prho, PE, h0, la, lap, lam;
    double u1g, u2g, u3g, u4g, rhog, ug, vg, pg;
    int i,j, k, l, Lid, Rid;
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
    Ap=row*col;

    row(0,0)=lt/a;
    row(1,0)=(u*lt/a)+nx*lc;
    row(2,0)=(v*lt/a)+ny*lc;
    row(3,0)=(h0*lt/a)+ud*lc;
    col(0,0)=-ud;
    col(0,1)=nx;
    col(0,2)=ny;
    col(0,3)=0;
    Ap+=row*col;
    for(j=0; j<4; j++){
        Ap(j,j)+=la;
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
    Am=row*col;

    row(0,0)=lt/a;
    row(1,0)=(u*lt/a)+nx*lc;
    row(2,0)=(v*lt/a)+ny*lc;
    row(3,0)=(h0*lt/a)+ud*lc;
    col(0,0)=-ud;
    col(0,1)=nx;
    col(0,2)=ny;
    col(0,3)=0;
    Am+=row*col;
    for(j=0; j<4; j++){
        Am(j,j)+=la;
    }
    Ap=Surface[surf].l*Ap;
    Am=Surface[surf].l*Am;
    x=Ap*Cell[Lid].U;
    if(Surface[surf].label==0){
        x+=Am*Cell[Rid].U;
        Cell[Lid].M=Cell[Lid].M+Ap;
        Cell[Lid].Res=Cell[Lid].Res-x;

        
        Cell[Rid].M=Cell[Rid].M-Am;
        Cell[Rid].Res=Cell[Rid].Res+x;
        for(i=0; i<4; i++){
            for(j=0; j<4; j++){
                A.insert(4*Lid+i,4*Rid+j)=Am(i,j);
                A.insert(4*Rid+i,4*Lid+j)=-Ap(i,j);
            }
        }
    }else{
        x+=Am*Ug;      
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
        Cell[Lid].M=Cell[Lid].M+Ap+Am*E;
        Cell[Lid].Res=Cell[Lid].Res-x;

    }

}

void Calculate_Flux(){
    int j;
    A.setZero();
    for(j=0; j<Ns; j++){
        Interface_Calculation(j);
    }
    printf("Finish Flux calculation\n");
}

void Update_U(){
    int i,j,k;
    double residual=1e5;
    

    //gather matrix
    Gather_Matrix();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver setup failed!" << std::endl;
    // return -1;
    }
    xi = solver.solve(bi);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        exit(1);
    }
    double temp=0;
    double max=0;
    for(i=0; i<Nc; i++){
        for(j=0; j<4; j++){
            Cell[i].U(j,0)=Cell[i].U(j,0)+xi(4*i+j);
        }
        temp+=Cell[i].Res(0)*Cell[i].Res(0)/Cell[i].vol*Cell[i].vol;
        Cell[i].rho=Cell[i].U(0,0);
        Cell[i].u=Cell[i].U(1,0)/Cell[i].rho;
        Cell[i].v=Cell[i].U(2,0)/Cell[i].rho;
        Cell[i].T=(Cell[i].U(3,0)/Cell[i].rho-0.5*(Cell[i].u*Cell[i].u+Cell[i].v*Cell[i].v))/CV;
        Cell[i].P=Cell[i].rho*R*Cell[i].T;
        Cell[i].a=sqrt(gamma*R*Cell[i].T);
        Cell[i].Res.setZero();
        Cell[i].M=Cell[i].vol/dt*Cell[i].M.setIdentity();
        
        max=((fabs(Cell[i].u)+Cell[i].a)>max)?(fabs(Cell[i].u)+Cell[i].a):max;
    }
    dt=(double)(CFL)/(max*double(Nx));
    printf("max=%f, res=%f dt=%f\n", max, sqrt(temp), dt);
    printf("Finish Update U\n");
}


void Gather_Matrix(){
    int j,k,l,m,n;
    int surf;
    for(j=0; j<Nc; j++){    
        for(k=0; k<4; k++){
            bi(4*j+k) = Cell[j].Res(k);
            surf=Cell[j].slist[k];
            for(l=0; l<4; l++){
                A.insert(4*j+k,4*j+l)=Cell[j].M(k,l);
            }
        }
    }
    A.makeCompressed(); 
}

void Save_Result(){
    int i,j;
    FILE* in;
    in = fopen("result.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            fprintf(in, "%lf %lf\n",Cell[i+j*Nx].rho, Cell[i+j*Nx].T);
        }
    }
    fclose(in);

    in = fopen("mesh.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            fprintf(in, "%lf %lf\n",Cell[i+j*Nx].x, Cell[i+j*Nx].y);
        }
    }
    fclose(in);

    in = fopen("u.txt","w");
    for(j=0; j<Ny; j++){
        for(i=0; i<Nx; i++){
            fprintf(in, "%lf %lf\n",Cell[i+j*Nx].u, Cell[i+j*Nx].v);
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