
/* to run, type
 
 gcc -Wall -o new_ss new_steady_state.c 
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define max_iterations 1000 


static int Rdim;
static double dr, V1, V2;
static double *radius;
static double *u_phi, *f, *U, *F;

FILE *output;
FILE *output1;
FILE *output2;
FILE *output3;



double RK(double value, double value2){
    double k1, k2, k3;
    
    k1 =value + (dr*value2);
    k2 = ((3.0/4.0)*value) + ((1.0/4.0)*k1) + ((dr/4.0)*value2);
    k3 = ((1.0/3.0)*value) + ((2.0/3.0)*k2) + ((2.0*dr/3.0)*value2);
    
    return k3;
    
}


double RK_with_update(double value, double value2, int position){
    double k1, k2, k3;
    double value_new;
    value_new = value2 - (value/(radius[position]*radius[position]));
    k1 =value + (dr*value_new);
    k2 = ((3.0/4.0)*value) + ((1.0/4.0)*k1) + ((dr/4.0)*(value2-(k1/(radius[position]))));
    k3 = ((1.0/3.0)*value) + ((2.0/3.0)*k2) + ((2.0*dr/3.0)*(value2-(k2/(radius[position]))));
    
    return k3;
    
}

void initialize(double r1){
    int i;
    radius[0] = r1;
    for(i=1; i<(Rdim+1); i++){
        radius[i] = radius[i-1] + dr;
    }

}

void find_s(double s){
    u_phi[0] = V1;
    f[0] = s;
    F[0] = 1.0;
    U[0] = 0.0;
    int i;
    for(i=1; i<(Rdim+1); i++){
        u_phi[i] = RK(u_phi[i-1], f[i-1]);
        f[i] = RK_with_update(f[i-1], (u_phi[i-1]/(radius[i-1]*radius[i-1])), (i-1));
        U[i] = RK(U[i-1], F[i-1]);
        F[i] = RK_with_update(F[i-1],(U[i-1]/(radius[i-1]*radius[i-1])), (i-1));
    }
    

}


int main(int argc, char **argv)
{
    
    Rdim = 100;
    double r1, r2;
    r1 = 1.0;
    r2 = 2.0;
    dr = (r2-r1)/Rdim;
    double omega;
    omega= 2.0;
    //V1 = r1*omega;
    //V2 = r2*omega;
    
    V1 = 10.0;
    V2 = 0.0;
    
    double s, s_old;
    
    s = 1.0;
    
    u_phi  = (double*) malloc((Rdim+1)*sizeof(double));
    f = (double*) malloc((Rdim+1)*sizeof(double));
    U = (double*) malloc((Rdim+1)*sizeof(double));
    F = (double*) malloc((Rdim+1)*sizeof(double));
    
    radius = (double*) malloc((Rdim+1)*sizeof(double));
    initialize(r1);
    
    find_s(s);
    printf("s = %f, u_phi[Rdim] = %f \n", s, u_phi[Rdim]);

    int i, g;
    for(i=0; i<max_iterations; i++){
        if(u_phi[Rdim-1]==V2){
            i=max_iterations;
        }else{
            s_old = s;
            s = s_old - ((u_phi[Rdim-1] - V2)/U[Rdim-1]);
            //printf("s_old = %f, s_new = %f, u_phi[RDIM] = %f, U[rdim] = %f \n", s_old, s, u_phi[Rdim], U[Rdim]);
            find_s(s);
        }
        
        //printf("s = %f, u_phi[Rdim] = %f \n", s, u_phi[Rdim]);
 
    }
    printf("s = %f, u_phi[Rdim] = %f \n", s, u_phi[Rdim-1]);
    
    output=fopen("Uphi.txt", "w");
    
    for(g = 0; g<(Rdim); g++)
    {
        fprintf(output, "%f \n", u_phi[g]);
        
    }




    return 0;
}

