
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


void open_file(){
    output2 = fopen("Uphi.txt", "rt");
    char line[80];
    int i=0;
    
    while(fgets(line, 80, output2) !=NULL){
        sscanf(line, "%lf",&d1uphi[i]);
        
        //printf("%f, i=%d \n", d1uphi[i], i);
        i+=1;
    }
    fclose(output2);
    
}


void fill_source(){
    int i;
    for(i=0; i<(P_size+2); i++){
        radius[i] = r1 + ((i-1)*grid_spacing[0]);
        if(i==0){
            source[i] = (2.0*radius[i]*Re*Wi*Wi);
        }else if(i==(P_size+1)){
            //source[i] = radius[i]*Re*Wo*Wo;
            source[i] = 0.0;
        }/*else if(i==1){ //this could be sketchy, not quite sure if this is the right thing to do.
            source[i] = Re*(d1uphi[i-1]/(radius));
        }else if(i==P_size){
            source[i] = 2.0*Re*(Wo*d1uphi[i-1]/(radius));
        }*/else{
            source[i] = Re*d1uphi[i-1]*d1uphi[i-1]/radius[i];
        }
        printf("source[%d] = %f, radius = %f \n ", i, source[i], radius[i]);

    }
    
}

void finite_difference(){
    int i;
    for(i=0; i<P_size; i++){
        //double radius = r1 + (i*grid_spacing[0]);
        double value = radius[i+1]*(fabs(source[i+2] - source[i]))/(grid_spacing[0]*Re);
        uphi_new[i] = sqrt(value);
        printf("value = %f, uphinew[%d] = %f \n", value, i, uphi_new[i]);
    }

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
    const int M = 102;
    const int N = 102;
    P_size = 100;
    int i, j;
    
    d1uphi = (double*) malloc(P_size*sizeof(double));
    uphi_new = (double*) malloc(P_size*sizeof(double));
    radius = (double*) malloc((P_size+2)*sizeof(double));


    source = (double*) malloc((P_size+2)*sizeof(double));
    
    open_file();
    Re = 1.0;
    Wi = 10.0;
    Wo = 0.0;
    r1 = 1.0;
    r2 = 2.0;
    
    
    grid_spacing[0] = (r2-r1)/P_size;
    grid_spacing[1] = (2.0*Pi)/P_size;
    
    fill_source();
    

    //printf("here \n");

    double *matrix2 = (double*) malloc(M*N*sizeof(double));
    // Declare an MxN matrix which can hold up to three band-diagonals.
    struct cs_sparse *triplet = cs_spalloc(M, N, 3*N, 1, 1);
    
    //cs_entry(triplet, 0, 0, 1.0); // Fill the corners with the value 1
    //cs_entry(triplet, M-1, N-1, 1.0);
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
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
    for (i=0; i<M; i++) {
        printf("%d \n", i);
        double a, b, c, d, e;
        //double radius = r1 + ((i-1)*grid_spacing[0]);
        //a = grid_spacing[1]*grid_spacing[1]*radius*(radius - grid_spacing[0]);
        //b = grid_spacing[0]*grid_spacing[0];
        //c = -2.0*((grid_spacing[1]*grid_spacing[1]*radius*radius) + (grid_spacing[0]*grid_spacing[0]));
        //d = b;
        //e = grid_spacing[1]*grid_spacing[1]*radius*(radius + grid_spacing[0]);
        //a = (1.0/(grid_spacing[0]*grid_spacing[0])) - (1.0/(radius*grid_spacing[0]));
        //b = 1.0/(radius*radius*grid_spacing[1]*grid_spacing[1]);
        //c = -2.0*((1.0/(grid_spacing[0]*grid_spacing[0])) + (1.0/(radius*radius*grid_spacing[1]*grid_spacing[1])));
        //c = -2.0*(1.0/(grid_spacing[0]*grid_spacing[0]));
        //d = b;
        //e = (1.0/(grid_spacing[0]*grid_spacing[0])) +(1.0/(radius*grid_spacing[0]));
        a = 1.0/grid_spacing[0];
        e = -a;
        
        /*if(i<=(P_size)){
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

        }*/
        
        
        if(i==0){
            //cs_entry(triplet, i, i, a);
            cs_entry(triplet, i, i+1, e+a);
        }else if(i==(M-1)){
            cs_entry(triplet, i, i-1, a+e);
            //cs_entry(triplet, i, i, e);
        }else{
            cs_entry(triplet, i, i-1, a);
            //cs_entry(triplet, i, i, c);
            cs_entry(triplet, i, i+1, e);
        }

        //cs_entry(triplet, i, i, c);
        //matrix2[(i*N) + i] = c;
    }
    
    // Declare and initialize the array of coefficients, 'b'.
    /*double *b = (double*) malloc(N*sizeof(double));
    double inner = 15.0;
    double outer =  -15.0;
    b[0] = inner;
    b[N-1] = outer;
    for (i=1; i<N-1; ++i) {
        //double radius = r1 + ((i-1)*grid_spacing[0]);
        //double coeffic = grid_spacing[0]*grid_spacing[0]*grid_spacing[1]*grid_spacing[1]*radius*radius;
        
        if((i-1)%P_size==0){
            b[i] = inner;
        }else if(i%P_size==0){
            b[i] = outer;
        }else{
            b[i] = 0.0;

        }
        
    }*/
    
    // Convert the triplet matrix into compressed column form, and use LU
    // decomposition to solve the system. Note that the vector 'b' of coefficients
    // overwritten to contain the solution vector 'x'.
    struct cs_sparse *matrix = cs_compress(triplet);
    cs_lusol(0, matrix, source, 1e-12);
    // Print the solution vector.
    output=fopen("trial3.txt", "w");

    for (i=0; i<(P_size+2); i++) {
        /*for(j=0; j<P_size+1; ++j){
            position = (i*P_size) + j;
            //printf(" m[%d] = %f", position, matrix2[position] );
            
            if(j==50){fprintf(output, "%+5.4e \n", b[i]);}

        }*/
        //printf("\n");
        printf("source[%d] = %+5.4e\n", i, source[i]);
        fprintf(output, "%+5.4e \n", source[i]);
    }
    
    for (i=0; i<N; ++i) {
        //printf("b[%d] = %+5.4e\n", i, b[i]);
    }
    
    
    //now want to finite difference the source, which is now the pressure and see if we get u back.
    //dP/dr = Re*Uphi^2/r

    finite_difference();
    
    
    output3=fopen("uphi_new.txt", "w");
    
    for (i=0; i<(P_size); i++) {
        //printf("\n");
        //printf("source[%d] = %+5.4e\n", i, source[i]);
        fprintf(output3, "%+5.4e \n", uphi_new[i]);
    }
    
    // Clean up memory usage.
    cs_spfree(triplet);
    cs_spfree(matrix);
    free(source);
    
    return 0;
}
