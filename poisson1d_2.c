
/* to run, type
 
gcc -Wall  -o poisson1d_2 poisson1d_2.c -I/$HOME/software/CSparse/Include -L/$HOME/software/CSparse/Lib -lcsparse
 
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cs.h"
#include "poisson1d_2.h"

FILE *output;
FILE *output1;
FILE *output2;
FILE *output3;

#define Pi 3.14159265
static double *uphi_new, *source, *radius;
static double Wi, Wo, r1, r2;
static double grid_spacing[2];
static int P_size;
static double Pinner;
//static double *uphi_new, *source, *radius;



/*void open_file(){
    output2 = fopen("Uphi.txt", "rt");
    output1 = fopen("Ur.txt", "rt");

    char line[80];
    int i=0;
    
    while(fgets(line, 80, output2) !=NULL){
        sscanf(line, "%lf",&d1uphi[i]);
        
        i+=1;
    }
    i=0;
    
    while(fgets(line, 80, output1) !=NULL){
        sscanf(line, "%lf",&d1ur[i]);
        
        i+=1;
    }
    fclose(output1);
    fclose(output2);
    
}*/

double deriv_r(int position, double *d1uphi){
    double value;
    value = (d1uphi[position] - d1uphi[position-P_size])/(grid_spacing[0]); 
    //printf("derivative in r of uphi = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}

double deriv_phi(int position, double *d1ur1){
    double value;

    double *d1ur  = &d1ur1[0];
    if(position%P_size==0){
        value = (d1ur[position+1] - d1ur[position+P_size - 1])/(grid_spacing[1]);
    }if((position+1)%P_size==0){
        //printf("here, position = %d, new position = %d \n", position, position-P_size +1);
        value = (d1ur[position-P_size + 1] - d1ur[position - 1])/(grid_spacing[1]);
    }else{
        value = (d1ur[position+1] - d1ur[position-1])/(grid_spacing[1]); 
    }
    
    //printf("derivative in phi of ur = %f, d1ur[%d] = %f \n", value, position, d1ur[position]);
    return value;
}

double delta_r_ur(int position, double *d1ur1){
    double value;
    double *d1ur  = &d1ur1[0];
    value = (d1ur[position] - d1ur[position-P_size])/(grid_spacing[0]); 
    //printf("derivative of ur in r = %E, d1ur[%d] = %E \n", value, position, d1ur[position]);
    return value;
}

void set_radius(){
    int i;
    for(i=0; i<(P_size+1); i++){
        radius[i] = r1 + ((i-1)*grid_spacing[0]);
    }
}

void fill_source(double *d1uphi, double *d1ur){
    int i, position;
    for(i=0; i<(P_size*(P_size+1)); i++){
        position = i/P_size;
        if(i<P_size){
            source[i] = Pinner;
        }else if(i>=P_size && i<(2*P_size)){
            source[i] = d1uphi[i-P_size]*d1uphi[i-P_size]/radius[position];
        }else{
            source[i] = 2.0*((deriv_r((i-P_size), d1uphi)*((d1uphi[i-P_size]/radius[position]) - (deriv_phi((i-P_size), d1ur)/radius[position]))) - (delta_r_ur((i-P_size), d1ur)*delta_r_ur((i-P_size), d1ur)));
        }
        //printf("source[%d] = %16.12e, radius = %f \n ", i, source[i], radius[position]);

    }
    
}

void finite_difference(){
    int i;
    for(i=0; i<(P_size); i++){
        double value = radius[i+1]*((source[((i+1)*P_size)] - source[i*P_size]))/(grid_spacing[0]);
        uphi_new[i] = value;
        //printf("value = %f, uphinew[%d] = %f \n", value, i, uphi_new[i]);
    }

}



void sparse(){
    
    const int M = P_size*(P_size+1);
    const int N = P_size*(P_size+1);
    int i;

    struct cs_sparse *triplet = cs_spalloc(M, N, 100*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double radius1 = radius[i/P_size];
        //printf("radius = %f, i = %d \n", radius1, i );
        double a, b, c, d, e, f, g;
        double deltar2 = 1.0/(grid_spacing[0]*grid_spacing[0]);
        double deltar_r = 1.0/(radius1*grid_spacing[0]);     
        double deltaphi = 1.0/(radius1*radius1*grid_spacing[1]*grid_spacing[1]);

        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = deltaphi;
        d = deltar2 + deltar_r - (2.0*deltaphi);
        e = deltaphi;
        //printf("i = %d, a = %f, b = %f, c = %f, d = %f, e = %f, radius = %f \n",i, a, b, c, d, e, radius1);
        
        
        f = (1.0/grid_spacing[0]);
        g = -1.0*f;
        
        if(i<P_size){
            cs_entry(triplet, i, i, 1.0);

        }else if(i>=P_size && i<(2*P_size)){
            
            
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else{
            
            if(i%P_size==0){
                cs_entry(triplet, i, i+(P_size-1), c);
                cs_entry(triplet, i, i+1, e);
            }else if((i+1)%P_size==0){
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i-(P_size-1), e);
            }else{
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i+1, e);
            }


            cs_entry(triplet, i, i-(2*P_size), a);
            cs_entry(triplet, i, i-(P_size), b);
            cs_entry(triplet, i, i, d);
            //printf("source[i-(2*Psize)] = %f, a = %f\n",source[i-(2*P_size)], a);
            //printf("source[i-(Psize)] = %f, c = %f\n",source[i-(P_size)], c );
            //printf("source[%d] = %f, e  = %f\n",i, source[i], e);

        }        
    }
    
    
    // Convert the triplet matrix into compressed column form, and use LU
    // decomposition to solve the system. Note that the vector 'b' of coefficients
    // overwritten to contain the solution vector 'x'.
    struct cs_sparse *matrix = cs_compress(triplet);
    //cs_print(matrix, 0);
    cs_lusol(0, matrix, source, 1e-10);
    cs_spfree(triplet);
    cs_spfree(matrix);
    
    
    

}


void poisson_pressure(double *d1uphi2, double *d1ur2, double *pressure, int P_size2, double r12, double r22, double Wi1, double Wo1)
// -----------------------------------------------------------------------------
// This program uses the CSparse library to solve the matrix equation A x = b.
// The values used for the coefficient 'b' and the matrix 'A' are totally
// arbitrary in this example. 'A' is given a band-tridiagonal form, with
// additional entries in the upper right and lower right corners.
// -----------------------------------------------------------------------------
{
    //for a 100x100 array, but have 1 ghost cell  for bcs
    //have a 100x100 pressure aray
    P_size = P_size2;
    Wi = Wi1;
    Wo = Wo1;
    r1 = r12;
    r2 = r22;
    int i;


    /*Wi = 5.0;
    Wo = 10.0;
    r1 = 1.0;
    r2 = 2.0;
    

    
    d1uphi = (double*) malloc(P_size*P_size*sizeof(double));
    
    d1ur = (double*) malloc(P_size*P_size*sizeof(double));*/


    source = (double*) malloc(P_size*(P_size+1)*sizeof(double));
    
    Pinner = 1.0; 
    radius = (double*) malloc((P_size+1)*sizeof(double)); 
    uphi_new = (double*) malloc(P_size*sizeof(double));

    //open_file();
    
    
    
    grid_spacing[0] = (r2-r1)/P_size;
    grid_spacing[1] = (2.0*Pi)/P_size;
    set_radius();
    
    fill_source(d1uphi2, d1ur2);
    sparse();

    
    
    // Print the solution vector.
    output=fopen("pressure.txt", "w");
    
    for (i=0; i<(P_size*(P_size+1)); i++) {

        //printf("source[%d] = %16.12e\n", i, source[i]);
        //pressure[i] = source[i];
        
        if(i%P_size==0){
            
            fprintf(output, "%f \n", source[i]);
        }
        
    }
    
    for (i=0; i<(P_size*(P_size+1)); i++) {
        //printf("source[%d] = %+5.4e\n", i, source[i]);
        pressure[i] = source[i];
        
    }
    
    
    finite_difference();
    output3=fopen("uphi_new.txt", "w");
    
    for (i=0; i<(P_size); i++) {
        
        fprintf(output3, "%f \n", uphi_new[i]);
    }
    
    
    
    fclose(output3);
    fclose(output);

    // Clean up memory usage.
    free(uphi_new);
    free(radius);
    free(source);
    
    //return 0;
}


