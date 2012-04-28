
/* to run, type
 
gcc -Wall  -o poisson1d_2 poisson1d_2.c -I/$HOME/software/CSparse/Include -L/$HOME/software/CSparse/Lib -lcsparse
 
 
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
static double Wi, Wo, r1, r2;
static double grid_spacing[2];
static int P_size;
static double Pinner;


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

double delta_r(int position){
    double value;
    if(position==1){
        value = (d1uphi[position-1] - (radius[position-1]*Wi))/(grid_spacing[0]);
    }else if(position==(P_size + 1)){       
        value = ((radius[position]*Wo) - d1uphi[position-2])/(grid_spacing[0]); 
    }else{
        value = (d1uphi[position-1] - d1uphi[position-2])/(grid_spacing[0]); 
    }
    //printf("derivative = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}

void set_radius(){
    int i;
    for(i=0; i<(P_size+1); i++){
        radius[i] = r1 + ((i-1)*grid_spacing[0]);
    }
}

void fill_source(){
    int i;
    for(i=0; i<(P_size+1); i++){
        
        if(i==0){
            source[i] = Pinner;
        }else if(i==1){

            source[i] = d1uphi[i-1]*d1uphi[i-1]/radius[i];
        }else{
            source[i] = 2.0*d1uphi[i-1]*delta_r(i)/radius[i];
        }
        printf("source[%d] = %f, radius = %f \n ", i, source[i], radius[i]);

    }
    
}

void finite_difference(){
    int i;
    for(i=0; i<(P_size); i++){
        double value = radius[i+1]*((source[i+1] - source[i]))/(grid_spacing[0]);
        uphi_new[i] = value;
        //printf("value = %f, uphinew[%d] = %f \n", value, i, uphi_new[i]);
    }

}



void sparse(){
    
    const int M = 101;
    const int N = 101;
    int i;
    // Declare an MxN matrix which can hold up to three band-diagonals.
    struct cs_sparse *triplet = cs_spalloc(M, N, 3*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double a, c, e, d, b;
        double deltar2 = 1.0/(grid_spacing[0]*grid_spacing[0]);
        double deltar_r = 1.0/(radius[i]*grid_spacing[0]);

        a = deltar2;
        c = -((2.0*deltar2) + (deltar_r));
        e = deltar2 + deltar_r;
        
        
        if(i==0){
            cs_entry(triplet, i, i, 1.0);

        }else if(i==1){
            d = 1.0/grid_spacing[0];
            b = -1.0*d;
            
            cs_entry(triplet, i, i-1, b);
            cs_entry(triplet, i, i, d);
            
        }else{
            cs_entry(triplet, i, i-2, a);
            cs_entry(triplet, i, i-1, c);
            cs_entry(triplet, i, i, e);
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
    //for a 100x100 array, but have 1 ghost cell  for bcs
    
    P_size = 100;
    int i;

    Wi = 5.0;
    Wo = 5.0;
    r1 = 1.0;
    r2 = 2.0;
    Pinner = 1.0;

    
    d1uphi = (double*) malloc(P_size*sizeof(double));
    uphi_new = (double*) malloc(P_size*sizeof(double));
    radius = (double*) malloc((P_size+1)*sizeof(double));


    source = (double*) malloc((P_size+1)*sizeof(double));
    
    open_file();
    
    
    
    grid_spacing[0] = (r2-r1)/P_size;
    grid_spacing[1] = (2.0*Pi)/P_size;
    set_radius();
    
    fill_source();
    
    sparse();

    
    
    // Print the solution vector.
    output=fopen("pressure.txt", "w");
    
    for (i=0; i<(P_size+1); i++) {
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
