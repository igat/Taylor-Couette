
/* to run, type
 
gcc -Wall  -o poisson poisson.c -I/$HOME/software/CSparse/Include -L/$HOME/software/CSparse/Lib -lcsparse
 
 
 */


#include <stdio.h>
#include <stdlib.h>
#include "cs.h"

FILE *output;

#define Pi 3.14159265

int main()
// -----------------------------------------------------------------------------
// This program uses the CSparse library to solve the matrix equation A x = b.
// The values used for the coefficient 'b' and the matrix 'A' are totally
// arbitrary in this example. 'A' is given a band-tridiagonal form, with
// additional entries in the upper right and lower right corners.
// -----------------------------------------------------------------------------
{
    //for a 10x10 array, but have 2 ghost cells in it for bcs
    const int M = 10002;
    const int N = 10002;
    const int P_size = 100;
    int i, j;
    //printf("here \n");

    double *matrix2 = (double*) malloc(M*N*sizeof(double));
    // Declare an MxN matrix which can hold up to three band-diagonals.
    struct cs_sparse *triplet = cs_spalloc(M, N, 6*N, 1, 1);
    
    cs_entry(triplet, 0, 0, 1.0); // Fill the corners with the value 1
    cs_entry(triplet, M-1, N-1, 1.0);
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    double grid_spacing[2];
    grid_spacing[0] = 1.0/M;
    grid_spacing[1] = (2*Pi)/N;
    double r1 = 1.0;
    int position;
    for(i=0; i<M; i++){
        for(j=0; j<N; j++){
            position = (i*N) + j;
            matrix2[position] = 0.0;
        }
    }
    
    matrix2[0] = 1.0;
    matrix2[(M-1)*N + N-1] = 1.0;
    
    //printf("here");
    for (i=1; i<M-1; ++i) {
        //printf("%d \n", i);
        double a, b, c, d, e;
        double radius = r1 + ((i-1)*grid_spacing[0]);
        a = grid_spacing[1]*grid_spacing[1]*radius*(radius - grid_spacing[0]);
        b = grid_spacing[0]*grid_spacing[0];
        c = -2.0*((grid_spacing[1]*grid_spacing[1]*radius*radius) + (grid_spacing[0]*grid_spacing[0]));
        d = b;
        e = grid_spacing[1]*grid_spacing[1]*radius*(radius + grid_spacing[0]);
        if(i<=(P_size)){
            cs_entry(triplet,i, 0, a);
            matrix2[(i*N)] = a;
        }else{
            cs_entry(triplet,i, i-P_size,  a);
            matrix2[(i*N) + (i-P_size)] = a;
        }
        
        if(i>=(M-P_size - 1)){
            cs_entry(triplet,i, N-1, e);
            matrix2[(i*N) + N-1] = e;
        }else{
            cs_entry(triplet,i, i+P_size,  e);
            matrix2[(i*N) + (i + P_size)] = e ;
        }
        
        
        if((i-1)%P_size==0){
            cs_entry(triplet,i, i+(P_size-1),  b);
            matrix2[(i*N) + (i + P_size-1)] = b;
        }else{
            cs_entry(triplet,i, i-1, b);
            matrix2[(i*N) + (i-1)] = b;
        }
        
        if(i%P_size==0){
            cs_entry(triplet,i, i-(P_size - 1), d);
            matrix2[(i*N) + (i - P_size +1)] = d;
        }else{
            cs_entry(triplet,i, i+1, d);
            matrix2[(i*N) + (i +1)] = d ;

        }

        cs_entry(triplet, i, i, c);
        matrix2[(i*N) + i] = c;
    }
    
    // Declare and initialize the array of coefficients, 'b'.
    double *b = (double*) malloc(N*sizeof(double));
    double inner = 0.5;
    double outer =  1.5;
    b[0] = inner;
    b[N-1] = outer;
    for (i=1; i<N-1; ++i) {
        if((i-1)%P_size==0){
            b[i] = inner;
        }else if(i%P_size==0){
            b[i] = outer;
        }else{
            b[i] = 0.0;

        }
        
    }
    
    // Convert the triplet matrix into compressed column form, and use LU
    // decomposition to solve the system. Note that the vector 'b' of coefficients
    // overwritten to contain the solution vector 'x'.
    struct cs_sparse *matrix = cs_compress(triplet);
    cs_lusol(0, matrix, b, 1e-12);
    // Print the solution vector.
    output=fopen("trial3.txt", "w");

    for (i=0; i<P_size+1; ++i) {
        for(j=0; j<P_size+1; ++j){
            position = (i*P_size) + j;
            //printf(" m[%d] = %f", position, matrix2[position] );
            
            if(j==50){fprintf(output, "%+5.4e \n", b[i]);}

        }
        //printf("\n");
        //printf("b[%d] = %+5.4e\n", i, b[i]);
    }
    
    for (i=0; i<N; ++i) {
        //printf("b[%d] = %+5.4e\n", i, b[i]);
    }

    
    
    
    
    // Clean up memory usage.
    cs_spfree(triplet);
    cs_spfree(matrix);
    free(b);
    
    return 0;
}
