
/* to run, type
 
gcc -Wall  -o poisson1d poisson1d.c -I/$HOME/software/CSparse/Include -L/$HOME/software/CSparse/Lib -lcsparse
 
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cs.h"

FILE *output;
FILE *output2;
FILE *output3;

#define Pi 3.14159265
static double *d1uphi, *uphi_new, *source, *radius;
static double Wi, Wo, Re, r1, r2;
static double grid_spacing[2];
static int P_size;
static double Pinner, Pouter;


void open_file(){
    output2 = fopen("Uphi.txt", "rt");
    char line[80];
    int i=0;
    
    while(fgets(line, 80, output2) !=NULL){
        sscanf(line, "%lf",&d1uphi[i]);
        
        i+=1;
    }
    fclose(output2);
    
}


void fill_source(){
    int i;
    for(i=0; i<(P_size+2); i++){
        //radius[i] = r1 + ((i-1)*grid_spacing[0]);
        if(i==0){
            //source[i] = (2.0*radius[i]*Re*Wi*Wi);
            source[i] = Re*radius[i+1]*Wi*Wi;
        }else if(i==(P_size+1)){
            //source[i] = radius[i]*Re*Wo*Wo;
            //source[i] = 0.0;
            //source[i] = Pouter;
            source[i] = Re*radius[i-1]*Wo*Wo;
        }else{
            source[i] = Re*d1uphi[i-1]*d1uphi[i-1]/radius[i];
        }
        printf("source[%d] = %f, radius = %f \n ", i, source[i], radius[i]);

    }
    
}

void finite_difference(){
    int i;
    for(i=0; i<P_size; i++){
        double value = radius[i+1]*(fabs(source[i+2] - source[i]))/(grid_spacing[0]*Re);
        uphi_new[i] = sqrt(value);
        //printf("value = %f, uphinew[%d] = %f \n", value, i, uphi_new[i]);
    }

}



void sparse(){
    
    const int M = 102;
    const int N = 102;
    int i;
    // Declare an MxN matrix which can hold up to three band-diagonals.
    struct cs_sparse *triplet = cs_spalloc(M, N, 3*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double a, e;
        
        a = 1.0/grid_spacing[0];
        e = -a;
        
        
        if(i==0){
            cs_entry(triplet, i, i, e);
            cs_entry(triplet, i, i+2, a);
        }else if(i==(M-1)){
            cs_entry(triplet, i, i, a);
            cs_entry(triplet, i, i-2, e);

        }else{
            cs_entry(triplet, i, i-1, a);
            cs_entry(triplet, i, i+1, e);
        }
        
        
    }
    
    
    // Convert the triplet matrix into compressed column form, and use LU
    // decomposition to solve the system. Note that the vector 'b' of coefficients
    // overwritten to contain the solution vector 'x'.
    struct cs_sparse *matrix = cs_compress(triplet);
    cs_lusol(0, matrix, source, 1e-12);
    cs_spfree(triplet);
    cs_spfree(matrix);

}


int main()
// -----------------------------------------------------------------------------
// This program uses the CSparse library to solve the matrix equation A x = b.
// The values used for the coefficient 'b' and the matrix 'A' are totally
// arbitrary in this example. 'A' is given a band-tridiagonal form, with
// additional entries in the upper right and lower right corners.
// -----------------------------------------------------------------------------
{
    //for a 100x100 array, but have 2 ghost cells in it for bcs
    
    P_size = 100;
    int i;
    //Pinner = 1.0;
    //Pouter = -10.337;
    
    Pinner = -10.0;
    Pouter = 10.0;
    
    d1uphi = (double*) malloc(P_size*sizeof(double));
    uphi_new = (double*) malloc(P_size*sizeof(double));
    radius = (double*) malloc((P_size+2)*sizeof(double));


    source = (double*) malloc((P_size+2)*sizeof(double));
    
    open_file();
    Re = 1.0;
    Wi = 5.0;
    Wo = 5.0;
    r1 = 1.0;
    r2 = 2.0;
    
    
    grid_spacing[0] = (r2-r1)/P_size;
    grid_spacing[1] = (2.0*Pi)/P_size;
    
    fill_source();
    
    sparse();
    
    /*if(Wi>Wo){
        if(source[0]>source[1] || source[P_size+1]>source[P_size]){
            if(source[0]>source[1]){
                Pinner = source[1];
            }
            
            if(source[P_size+1]>source[P_size]){
                Pouter = source[P_size];
            }
            fill_source();
            sparse();
            
        }
    }else{
        printf("here, Wi = %f, Wo = %f \n", Wi, Wo);
        if(source[0]<source[1] || source[P_size+1]<source[P_size]){
            printf("source[0] = %f, source[1] = %f, source[P+1] = %f, souce[P] = %f \n", source[0], source[1], source[P_size+1] , source[P_size]);
            if(source[0]<source[1]){
                Pinner = source[1];
            }
            
            if(source[P_size+1]<source[P_size]){
                Pouter = source[P_size];
            }
            fill_source();
            sparse();
        }
    }*/
    


    // Print the solution vector.
    output=fopen("pressure.txt", "w");

    for (i=0; i<(P_size+2); i++) {
        printf("source[%d] = %+5.4e\n", i, source[i]);
        fprintf(output, "%+5.4e \n", source[i]);
    }
    
    finite_difference();
    
    
    output3=fopen("uphi_new.txt", "w");
    
    for (i=0; i<(P_size); i++) {

        fprintf(output3, "%+5.4e \n", uphi_new[i]);
    }
    
    // Clean up memory usage.

    free(source);
    
    return 0;
}
