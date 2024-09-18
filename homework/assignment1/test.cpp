//Problem 5.5 (Elliptic equation solution around thin airfoil)
#include<iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<vector>
#include<math.h>

using namespace std;

//f(kappa)
long double f_kappa (long double dy_min, long double k, int jL, long double D)
{
    long double ktop = k/(static_cast<long double>(jL-1));
    return dy_min - D*( (exp(ktop)-1) / (exp(k)-1) );
}

// Derivative of f(kappa) with respect to k
long double f_dash_kappa(long double k, int jL, long double D)
{
    long double N = static_cast<long double>(jL-1);
    return  D*( N*(exp(k/N) - 1)*exp(k) - (exp(k)-1)*exp(k/N) ) / ( N*(exp(k)-1)*(exp(k)-1) );
}

//compute  kappa param
long double get_kappa(long double dy_min, long double k_init, int jL, long double D)
{
    long double kappa = k_init;
    for (int it = 0; it<10; it++)
    {
        kappa = kappa - f_kappa(dy_min, kappa, jL, D)/f_dash_kappa(kappa, jL, D);
         std::cout<<kappa<<"\n";
    }
    return kappa;
}

//get residual given a vector of phi and A
long double get_residual ( vector<vector<long double> > phi_new,
                           vector<vector<long double> > x,
                           vector<vector<long double> > y, 
                           long double A)
{
    size_t N_t = phi_new[0].size();
    long double residual =0.0;
    for (int i=1; i<N_t-1; i++)
    {
        for (int j=1; j<N_t-1; j++)
        {
            long double diff1 = (phi_new[i+1][j] - phi_new[i][j])/(x[i+1][j]-x[i][j]);
            long double diff2 = (phi_new[i][j] - phi_new[i-1][j])/(x[i][j]-x[i-1][j]);
            long double A_dphi2dx2 = 2*A*(diff1-diff2)/(x[i+1][j]-x[i-1][j]);

            diff1 = (phi_new[i][j+1] - phi_new[i][j])/(y[i][j+1]-y[i][j]);
            diff2 = (phi_new[i][j] - phi_new[i][j-1])/(y[i][j]-y[i][j-1]);
            long double dphi2dy2 = (diff1-diff2)/(y[i][j+1]-y[i][j-1]);
            // cout<<abs(A_dphi2dx2+dphi2dy2)<<"\n";
            residual = residual> abs(A_dphi2dx2+dphi2dy2)? residual : abs(A_dphi2dx2+dphi2dy2);
        }
    }
    return residual;
}

vector<long double> calculate_cp(
    const vector<vector<long double> >& phi_new,
    const vector<vector<long double> >& x,
    long double dy_min,
    long double p_inf,
    long double g,
    long double M_inf,
    long double V_inf,
    long double rho_inf)
{
    size_t N_t = phi_new.size();
    vector<long double> cp(N_t, 0.0);
    long double u, v;

    for (size_t i = 0; i < N_t; i++)
    {
        // fill u
        if (i == 0)
        {
            u = (phi_new[i+1][0] - phi_new[i][0]) / (x[i+1][0] - x[i][0]);
        }
        else if (i == N_t - 1)
        {
            u = (phi_new[i][0] - phi_new[i-1][0]) / (x[i][0] - x[i-1][0]);
        }
        else
        {
            u = (phi_new[i+1][0] - phi_new[i-1][0]) / (x[i+1][0] - x[i-1][0]);
        }

        // fill v
        v = (phi_new[i][1] - phi_new[i][0]) / dy_min;
        long double p = p_inf * pow((1 - 0.5 * (g - 1) * M_inf * M_inf * ((u * u + v * v) / (V_inf * V_inf) - 1)), g / (g - 1));
        long double q = 0.5 * rho_inf * V_inf * V_inf;
        cp[i] = (p - p_inf) / q;
    }

    return cp;
}

//main program
int main()
{
    //list of constants
    const int N_c = 21; // Number of points across the chord
    const int N_o = 15; // Number of points outside the chord
    const int N_t = 51; // Total number of points along each direction
    const int method = 2; // which method to solve
    const int nsteps = 1000; // number of iterations
    const long double c = 1.0; // Chord length
    const long double L_x = 50*c; // Length along x outside chord until freestream
    const long double L_y = 50*c;  // Length along y 

    //more constants derived for the problem for grid stretching along x
    const int i_S = 0; // starting index at freestream
    const int i_L = 15; // index for the leading edge
    const int i_T = i_L + N_c-1; // index for the trailing edge
    const long double x_S = -50; // starting location for the grid
    const long double y_S = 0;
    const long double x_L = 0; // leading edge x-coord
    const long double x_T = 1; // trailing edge x-coord
    const long double x_E = 51; // freestream end x-coord
    const long double dx_o = L_x/(N_o-1); // dx outside the chord in freestream for equispaced mesh
    const long double dx_c = c/(N_c-1); // dx inside the chord
    const long double dy = L_y/(N_t-1); // mesh spacing along y 
    const long double th = 0.06; // thickness of airfoil
    const long double dy_min = th/10; // first mesh point spacing along y 
    const long double kappa_y = get_kappa(dy_min, 1.0, N_t, L_y); // mesh stretching params along y 
    const long double kappa_x_L = get_kappa(dx_c, 1.0, N_o-1, L_x-dx_c); // mesh stretching params along x
    const long double kappa_x_T = get_kappa(dx_c, 1.0, N_o-1, L_x-dx_c); // mesh stretching params along x

    // physical params start here
    const long double g = 1.4; // gamma = cp/cv
    const long double p_inf=1; // freestream pressure
    const long double rho_inf=1; // freestream density
    const long double M_inf = 0.5; //freestream mach number
    const long double A = 1-M_inf*M_inf; // constant for x-der
    const long double a_inf = sqrt(g*p_inf/rho_inf); // freestream speed of sound 
    const long double V_inf = M_inf*a_inf; // freestream velocity

    // Declare a 2D vector of size 51x51
    vector<vector<long double> > x(N_t, vector<long double>(N_t, 0.0)); // x-coords
    vector<vector<long double> > y(N_t, vector<long double>(N_t, 0.0)); // y-coords
    vector<vector<long double> > phi_old(N_t, vector<long double>(N_t, 0.0)); // \phi^n
    vector<vector<long double> > phi_new(N_t, vector<long double>(N_t, 0.0)); // \phi^{n+1}

    //output file
    ofstream gridFile("initial_grid.txt");
    ofstream icFile("initial_condition.txt");
    ofstream resFile("residuals.txt");
    ofstream cpFile("pressurecoeff.txt");
    ofstream jobFile("jobinfo.txt");
    jobFile<<method<<"\n";
    jobFile.close();

    //setup 
    for (size_t j=0; j<N_t; j++)
    {
        for (size_t i=i_S; i<N_t; i++)
        {
            if(i<i_L) // before leading edge
            {
                //equispaced
                // x[i][j] = x_S + (i)*dx_o;

                //stretched
                long double ktop = kappa_x_L * (static_cast<long double>(i_L-i-1) / static_cast<long double>(N_o-1));
                x[i][j] = x_L-dx_c - (L_x-dx_c)*( (exp(ktop)-1) / (exp(kappa_x_L)-1) );
            }
            else if (i>=i_L && i<=i_T) // inside the chord
            {
                x[i][j] = x_L + (i-i_L)*dx_c; //equispaced
            }
            else // outside trailing edge
            {
                // equispaced 
                // x[i][j] = x_T + (i-i_T)*dx_o; 

                //stretched
                long double ktop = kappa_x_T * (static_cast<long double>(i-i_T-1) / static_cast<long double>(N_o - 1));
                x[i][j] = x_T+dx_c + (L_x-dx_c)*( (exp(ktop)-1) / (exp(kappa_x_T)-1) );
            }

            //now compute stretched y-coords
            long double k_top = kappa_y * (static_cast<long double>(j) / static_cast<long double>(N_t -1));
            y[i][j] = y_S + L_y*( (exp(k_top)-1) / (exp(kappa_y)-1) ); 
        }
    }

    //Write out grid file
    for (size_t i=0; i<N_t; i++)
    {
        for (size_t j=0; j<N_t; j++)
        {
            gridFile<<x[i][j]<<",\t"<<y[i][j]<<"\n";
        }
    }
    gridFile.close();

    //The actual solution begins from here. 
    //We have the grids. Next, set the initial condition
    for (size_t i=0; i<N_t; i++)
    {
        for (size_t j=0; j<N_t; j++)
        {
            phi_old[i][j] = V_inf*x[i][j];
            icFile<<x[i][j]<<",\t"<<y[i][j]<<",\t"<<phi_old[i][j]<<"\n";
        }
    }

    //start the outer iteration loop
    for (int n=0; n<nsteps; n++)
    {
        // left, right and top boundaries depend only on phi_old so compute them first
        //left and right boundaries 
        for (int j=0; j<N_t; j++)
        {
            if (method==1)
            {
                phi_new[0][j] = V_inf*x[0][j];
                phi_new[N_t-1][j] = V_inf*x[N_t-1][j];
            }
            else if (method==2)
            {
                phi_old[0][j] = V_inf*x[0][j];
                phi_old[N_t-1][j] = V_inf*x[N_t-1][j];
            }
        }

        // top boundary
        for (int i=0;i<N_t; i++)
        {
            if (method==1)
            {
                phi_new[i][N_t-1] = V_inf*x[i][N_t-1];
            }
            else if (method==2)
            {
                phi_old[i][N_t-1] = V_inf*x[i][N_t-1];
            }
        }
        
        //interior points
        for (int i=1; i<N_t-1; i++)
        {
            for (int j=1; j<N_t-1; j++)
            {
                long double alpha(0.0), beta(0.0), gamma(0.0), delta(0.0);
                if (method==1) // point-jacobi
                {
                    alpha = 2.0*A/(x[i+1][j]-x[i][j])/(x[i+1][j]-x[i-1][j]);
                    beta  = 2.0*A/(x[i][j]-x[i-1][j])/(x[i+1][j]-x[i-1][j]);
                    gamma = 2.0/(y[i][j+1]-y[i][j])/(y[i][j+1]-y[i][j-1]);
                    delta = 2.0/(y[i][j]-y[i][j-1])/(y[i][j+1]-y[i][j-1]);
                    phi_new[i][j] = (1/(alpha+beta+gamma+delta))*( alpha * phi_old[i+1][j] 
                                                                 + beta  * phi_old[i-1][j] 
                                                                 + gamma * phi_old[i][j+1] 
                                                                 + delta * phi_old[i][j-1] );
                }
                else if (method==2) //point-Gauss-Siedel
                {
                    alpha = 2.0*A/(x[i+1][j]-x[i][j])/(x[i+1][j]-x[i-1][j]);
                    beta  = 2.0*A/(x[i][j]-x[i-1][j])/(x[i+1][j]-x[i-1][j]);
                    gamma = 2.0/(y[i][j+1]-y[i][j])/(y[i][j+1]-y[i][j-1]);
                    delta = 2.0/(y[i][j]-y[i][j-1])/(y[i][j+1]-y[i][j-1]);
                    phi_old[i][j] = (1/(alpha+beta+gamma+delta))*( alpha * phi_old[i+1][j] 
                                                                 + beta  * phi_old[i-1][j] 
                                                                 + gamma * phi_old[i][j+1] 
                                                                 + delta * phi_old[i][j-1] );

                }
            }
        }

         //bottom boundary needs latest data, fill it after interior points
        for (int i=0; i<N_t; i++)
        {
            
            //outside the airfoil surface: zero gradient in y
            if(i<=i_L || i>=i_T)
            {
                if (method==1)
                {
                    phi_new[i][0] = phi_new[i][1];
                }
                else if (method==2)
                {
                    phi_old[i][0] = phi_old[i][1];
                }
            }
            //on the airfoil surface, impose the slope-based B.C. 
            else
            {
                //get slope
                //double dydx = (th/2/sqrt(x[i][0]*(c-x[i][0])))*(1-2*(x[i][0]/c)); //circular
                double dydx = -1.0*(x[i][0] -  c/2)/sqrt(( pow(c*c/4/th+th/4,2) -pow(x[i][0]-c/2,2))); //elliptic
                if (method==1)
                {phi_new[i][0] = phi_new[i][1]-V_inf*dy_min*dydx;}
                else if (method==2)
                {phi_old[i][0] = phi_old[i][1]-V_inf*dy_min*dydx;}
            }
        }

        //compute residuals (MSE residual between iters)
        // long double sum[nsteps] ={0.0};
        // for (size_t i=0;i<N_t;i++)
        // {
        //     for (size_t j=0; j<N_t; j++)
        //     {
        //         sum[n] = sum[n] + abs(phi_old[i][j]-phi_new[i][j])*abs(phi_old[i][j]-phi_new[i][j]);
        //     }
        // }
        // sum[n] = sqrt(sum[n])/N_t;
        // cout<<sum[n]<<"\n";

        //equation residual
        long double residual[nsteps]={0.0};
        if (method==1)
        {
            residual[n] = get_residual(phi_new, x, y, A);
        }
        else if (method==2)
        {
            residual[n] = get_residual(phi_old, x, y, A);
        }
        

        //write out residual
        resFile<<n<<",\t"<<residual[n]<<"\n";

        //update phi_old to phi_new for point jacobi
        if (method==1)
        {
            for (size_t i=0; i<N_t; i++)
            {
                for (size_t j=0; j<N_t; j++)
                {
                    phi_old[i][j] = phi_new[i][j];
                }
            }
        }
    }

    // get pressure coefficient
    vector<long double> cp;
    if (method==1)
    {
        //Point Jacobi
        cp = calculate_cp(phi_new, x, dy_min, p_inf, g, M_inf, V_inf, rho_inf);
    }
    else if (method==2)
    {
        //Point Gauss Siedel
        cp = calculate_cp(phi_old, x, dy_min, p_inf, g, M_inf, V_inf, rho_inf);
    }

    // Write cp to file
    for (size_t i = 0; i < N_t; i++)
    {
        cpFile << x[i][0] << ",\t" << cp[i] << "\n";
    }
    //close output files
    resFile.close();
    cpFile.close();

    return 0;
}
