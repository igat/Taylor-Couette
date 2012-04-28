
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

double delta_r(int position){
    double value;
    if(position==1){     //reflecting boundary
        value = (d1uphi[position-1] - (radius[position-1]*Wi))/(grid_spacing[0]);
        //value = Wi;
        //value = (2.0*d1uphi[position])/(grid_spacing[0]);
    }else if(position==(P_size)){        //reflecting boundary everywhere
        value = ((radius[position+1]*Wo) - d1uphi[position-2])/(grid_spacing[0]); 
        //value = Wo;
        //value = (-2.0*d1uphi[position-2])/(grid_spacing[0]); 
    }else{
        value = (d1uphi[position-1] - d1uphi[position-2])/(grid_spacing[0]); 
    }
    printf("derivative = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}
/*
double delta_r(int position){
    double value;
    if(position==1){     //reflecting boundary
        value = ((radius[position + 1]*d1uphi[position]) - (radius[position-1]*radius[position-1]*Wi))/(radius[position]*grid_spacing[0]);
        //value = Wi;
    }else if(position==(P_size)){        //reflecting boundary everywhere
        value = ((radius[position+1]*radius[position+1]*Wo) - (radius[position-1]*d1uphi[position-2]))/(radius[position]*grid_spacing[0]); 
        //value = Wo;
    }else{
        value = ((radius[position + 1]*d1uphi[position]) - (radius[position-1]*d1uphi[position-2]))/(radius[position]*grid_spacing[0]); 
    }
    printf("derivative = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}*/



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
            //source[i] = (2.0*radius[i]*Re*Wi*Wi);
            source[i] = Pinner;
            //source[i] = Re*radius[i+1]*Wi*Wi;
        }else if(i==(P_size+1)){
            //source[i] = radius[i]*Re*Wo*Wo;
            //source[i] = 0.0;
            source[i] = radius[i]*Wo*Wo;
        }else{
            source[i] = 2.0*d1uphi[i-1]*delta_r(i)/radius[i];
        }
        //source[i] = 2.0*Re*d1uphi[i]*delta_r(i)/radius[i];
        printf("source[%d] = %f, radius = %f \n ", i, source[i], radius[i]);

    }
    
}

void finite_difference(){
    int i;
    //uphi_new[0] = r1*Wi;
    //uphi_new[P_size-1] = Wo*(r2-grid_spacing[0]);
    for(i=0; i<(P_size); i++){
        double value = radius[i+1]*(fabs(source[i+1] - source[i]))/(grid_spacing[0]);
        uphi_new[i] = sqrt(value);
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
        a = (deltar2 - deltar_r);
        c = -2.0*deltar2;
        e = (deltar2 + deltar_r);
        //d = 1.0/grid_spacing[0];
        //b = -1.0*d;
        //d = radius[i+2]/(radius[i+1]*grid_spacing[0]);
        //b = -1.0*radius[i]/(radius[i+1]*grid_spacing[0]);
        
        if(i==0){
            //d = radius[i+2]/(radius[i+1]*grid_spacing[0]);
            //b = -1.0*radius[i]/(radius[i+1]*grid_spacing[0]);
            cs_entry(triplet, i, i, 1.0);
            //cs_entry(triplet, i, i+2, d);
            //cs_entry(triplet, i, i+1, c);
            //cs_entry(triplet, i, i, a);
            //cs_entry(triplet, i, i+2, e);
        }else if(i==(M)){
            d = 1.0/grid_spacing[0];
            b = -1.0*d;
            
            cs_entry(triplet, i, i-1, b);
            cs_entry(triplet, i, i, d);
            /*cs_entry(triplet, i, i-2, a);
            cs_entry(triplet, i, i, e);
            cs_entry(triplet, i, i-1, c);*/
            
        }else{
            cs_entry(triplet, i, i-1, a);
            cs_entry(triplet, i, i, c);
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
    Re = 1.0;
    Wi = 5.0;
    Wo = 5.0;
    r1 = 1.0;
    r2 = 2.0;
    Pinner = 1.0;
    //Pouter = 40.0;
    //Pinner = Re*r1*r1*Wi*Wi/2.0;
    //Pouter = ((Re*r2*r2*Wo*Wo)/2.0) + ;
    
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
