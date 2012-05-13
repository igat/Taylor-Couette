
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
static double *uphi_new, *source, *source_old, *source_center, *radius;
static double Wi, Wo, r1, r2, Re;
static double grid_spacing[2];
static int P_size;
static double Pinner;
//static double *uphi_new, *source, *radius;
#define TRUE  1
#define FALSE 0
typedef int boolean;
static double pertubs;

double deriv_r(int position, double *d1uphi){
    double value;
     double *w  = &d1uphi[0];
    
    if(position<=(P_size-1)){     //reflecting boundary
        value = Wi;
    }else if(position>=(((P_size-1)*P_size))){// && position<(((P_size)*P_size)-7)){        //reflecting boundary everywhere
        value = Wo;
        //value = ((Wo*r2) - w[position-P_size])/(2.0*grid_spacing[0]);
    }else{
        //value = (w[position + N[1]] - w[position-N[1]])/(grid_spacing[0]); 
        value = (w[position+P_size] - w[position-P_size])/(2.0*grid_spacing[0]); 
        //value = (w [position] - w [position-P_size])/(grid_spacing[0]); 
    }
    value = (w [position] - w [position-P_size])/(grid_spacing[0]);
    //printf("derivative in r of uphi = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}

double deriv_phi(int position, double *d1ur1){
    double value;
    double *w  = &d1ur1[0];
    /*if(position%P_size==0){
        //value = (w[position + 1] - w[position + P_size-1])/(2.0*grid_spacing[1]);
        value = (w[position] - w[position+ P_size-1])/(grid_spacing[1]);
    }else if(position%P_size==(P_size-1)){
        //value = (w[position-P_size+1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position] - w[position-1])/(grid_spacing[1]);
    }else{
        //value = (w[position + 1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position] - w[position-1])/(grid_spacing[1]);
    }
    if(position%P_size==(P_size-1)){
        //value = (w[position + 1] - w[position + P_size-1])/(2.0*grid_spacing[1]);
        value = (w[position-(P_size-1)] - w[position])/(grid_spacing[1]);
    }else if(position%P_size==(P_size-1)){
        //value = (w[position-P_size+1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position+1] - w[position])/(grid_spacing[1]);
    }else{
        //value = (w[position + 1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position+1] - w[position])/(grid_spacing[1]);
    }*/
    
    if(position%P_size==0){
        //value = (w[position + 1] - w[position + P_size-1])/(2.0*grid_spacing[1]);
        value = (w[position+1] - w[position+ P_size-1])/(2.0*grid_spacing[1]);
    }else if(position%P_size==(P_size-1)){
        //value = (w[position-P_size+1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position-(P_size-1)] - w[position-1])/(2.0*grid_spacing[1]);
    }else{
        //value = (w[position + 1] - w[position-1])/(2.0*grid_spacing[1]);
        value = (w[position+1] - w[position-1])/(2.0*grid_spacing[1]);
    }

    return value;
}


double delta_r_ur(int position, double *d1ur1){
    double value;
    double *w  = &d1ur1[0];
    if(position<=(P_size-1)){     //reflecting boundary
        value = (w[position + P_size])/(grid_spacing[0]);
        //value = (w [position])/(grid_spacing[0]);

        //value = (2.0*w[position + N[1]])/(grid_spacing[0]);
    }else if(position>=((P_size-1)*P_size)){        //reflecting boundary everywhere
        //value = (-w[position-P_size])/(2.0*grid_spacing[0]); 
        //value = (w [position] - w [position-P_size])/(grid_spacing[0]);

        value = (-w[position-P_size])/(grid_spacing[0]); 
        //value = 0.0;
    }else{
        //value = (w[position + N[1]] - w[position-N[1]])/(grid_spacing[0]); 
        value = (w[position+P_size] - w[position-P_size])/(2.0*grid_spacing[0]); 
        //value = (w [position] - w [position-P_size])/(grid_spacing[0]); 
    }
    value = (w [position] - w [position-P_size])/(grid_spacing[0]);

    
    //printf("derivative of ur in r = %E, d1ur[%d] = %E \n", value, position, d1ur[position]);
    return value;
}

void set_radius(){
    int i;
    for(i=0; i<(P_size+2); i++){
        radius[i] = r1 + ((i-1)*grid_spacing[0]);
    }
}
void deriv_phi_new(double *new_uphi, double *old_uphi){
    int i, j, position2;
    for(i=0; i<P_size; i++){
        for(j=0; j<P_size; j++){
            position2 = (i*P_size) + j;
            new_uphi[position2] = deriv_r(position2,old_uphi);
        }
    }
    
}


double deriv_phi2(double *vel, int position){
    double value;
    double *w  = &vel[0];
    /*if(position%P_size==0){
        value = (w[position] + w[position + P_size -2] - (2.0*w[position + P_size-1]))/(grid_spacing[1]*grid_spacing[1]);
        //printf("at boundary condition (0) for phi2. position = %d \n", position);
        //value = (w[position + 1] + w[position + P_size-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else if(position%P_size==1){
        value = (w[position] + w[position+P_size-2] - (2.0*w[position-1]))/(grid_spacing[1]*grid_spacing[1]);
    }else{
        //value = (w[position + 1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
        value = (w[position] + w[position-2] - (2.0*w[position-1]))/(grid_spacing[1]*grid_spacing[1]);
    }*/
    /*if(position%P_size==(P_size-1)){
        value = (w[position- (P_size-2)] + w[position] - (2.0*w[position - (P_size-1)]))/(grid_spacing[1]*grid_spacing[1]);
        //printf("at boundary condition (0) for phi2. position = %d \n", position);
        //value = (w[position + 1] + w[position + P_size-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else if(position%P_size==(P_size-2)){
        value = (w[position- (P_size-2)] + w[position] - (2.0*w[position+1]))/(grid_spacing[1]*grid_spacing[1]);
    }else{
        //value = (w[position + 1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
        value = (w[position+1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }*/
    
    if(position%P_size==0){
        value = (w[position+1] + w[position+(P_size-1)] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
        //printf("at boundary condition (0) for phi2. position = %d \n", position);
        //value = (w[position + 1] + w[position + P_size-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else if(position%P_size==(P_size-1)){
        value = (w[position- (P_size-1)] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else{
        //value = (w[position + 1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
        value = (w[position+1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }

    
    
    
    return value;
}

void fill_source(double *d1uphi, double *d1ur){
    int i, position;
    double *value_new = (double*) malloc(P_size*P_size*sizeof(double));
    deriv_phi_new(value_new, d1uphi);
    
    for(i=0; i<(P_size*(P_size+2)); i++){
        position = i/P_size;
        if(i<(P_size)){
            source[i] = Pinner;
            //printf("pinner = %f \n", Pinner);
        }else if(i>=(P_size) && i<(2*P_size)){
            //source[i] = d1uphi[i-P_size]*d1uphi[i-P_size]/radius[position];
            source[i] = r1*Wi*Wi;
        }else if(i>=(P_size*(P_size+1)) ){
            source[i] = r2*Wo*Wo;
        }else{
            source[i] = 2.0*((deriv_r((i-P_size), d1uphi)*((d1uphi[i-P_size]/radius[position]) - (deriv_phi((i-P_size), d1ur)/radius[position]))) - (delta_r_ur((i-P_size), d1ur)*delta_r_ur((i-P_size), d1ur)) - ((1.0/(Re*radius[position]*radius[position]))*(delta_r_ur((i-P_size), d1ur) + (deriv_phi((i-P_size), value_new)) - (deriv_phi2(d1ur, (i-P_size))/radius[position]))));
        }
        //source_old[i] = source[i];
        //source_center[i] = source[i];
        //printf("source[%d, %d, %d] = %f \n", i , position, i%P_size, source_old[i] );

    }
    free(value_new);
    
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
    
    const int M = (P_size*(P_size+2));
    const int N = (P_size*(P_size+2));
    int i;

    struct cs_sparse *triplet = cs_spalloc(M, N, 100*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double radius1 = radius[i/(P_size)];
        //printf("radius = %f, i = %d, radius[i-1] = %f \n", radius1, i, radius[(i/P_size)-1] );
        double a, b, c, d, e1, f, g;
        double deltar2 = 1.0/(grid_spacing[0]*grid_spacing[0]);
        double deltar_r = 1.0/(radius1*grid_spacing[0]);     
        double deltaphi = 1.0/(radius1*radius1*grid_spacing[1]*grid_spacing[1]);
/*
        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = deltaphi;
        d = deltar2 + deltar_r - (2.0*deltaphi);
        e1 = deltaphi;*/
        //printf("i = %d, a = %f, b = %f, c = %f, d = %f, e = %f, radius = %f \n",i, a, b, c, d, e, radius1);
        /*a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r;// + c ;
        e1 = deltaphi;*/
        /*
        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r + deltaphi ;
        e1 = deltaphi;*/

        /*
        a = deltar2*(radius[(i/P_size)-1]/radius1);
        b = ((-((radius[(i/P_size)-1] + radius1)/radius1)*deltar2) - (deltar_r));
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r + deltaphi;
        e1 = deltaphi;*/
        /*
        a = deltar2*(radius[(i/(P_size))-1]/radius1);
        b = ((-(radius[(i/(P_size))-1]/radius1)*deltar2) - (deltar2));
        c = -2.0*deltaphi;
        d = deltar2;// + deltaphi;
        e1 = deltaphi;*/
        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r + deltaphi;
        e1 = deltaphi;
        
        
        f = (1.0/(grid_spacing[0]));
        g = -1.0*f;
        if(i<(P_size)){
            cs_entry(triplet, i, i, 1.0);
            
        }else if(i>=(P_size) && i<(2*P_size)){
            
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else if(i>=((P_size)*(P_size+1))){
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else{

            if(i%P_size==(P_size-1)){
                cs_entry(triplet, i, i-(P_size-1),c);// c+e1);
                cs_entry(triplet, i, i-(P_size-2), e1);
            }else if(i%P_size==(P_size-2)){
                cs_entry(triplet, i, i+1,c);// c+e1);
                cs_entry(triplet, i, i-(P_size-2), e1);
            }else{
                cs_entry(triplet, i, i+1, c);
                cs_entry(triplet, i, i+2, e1);
            }

            
            /*
            if(i%P_size==0){
                cs_entry(triplet, i, i+(P_size-1), c);
                cs_entry(triplet, i, i+(P_size-2), e1);
            }else if(i%P_size==1){
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i+(P_size-2), e1);
            }else{
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i-2, e1);
            }*/
            /*
            if(i%P_size<=P_size/4){
                //cs_entry(triplet, i, i+P_size-1, e1);
                //cs_entry(triplet, i, i+1, e1);
                d = deltar2 + deltar_r + deltaphi;
                cs_entry(triplet, i, i+1, c);
                cs_entry(triplet, i, i+2, e1);
            }else if(i%P_size==(P_size-1)){
                //cs_entry(triplet, i, i-1, e1 );
                //cs_entry(triplet, i, i-P_size+1, e1);
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i-2, e1);
                d = deltar2 + deltar_r + deltaphi;

            }else{
                cs_entry(triplet, i, i-1, e1);
                cs_entry(triplet, i, i+1, e1);
                 d = deltar2 + deltar_r+c;

            }*/
           
            
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
    cs_lusol(0, matrix, source, 1e-15);
    cs_spfree(triplet);
    cs_spfree(matrix);
    
    
    

}
/*
void sparse2(){
    
    const int M = (P_size*(P_size+2));
    const int N = (P_size*(P_size+2));
    int i;
    
    struct cs_sparse *triplet = cs_spalloc(M, N, 100*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double radius1 = radius[i/(P_size)];
        //printf("radius = %f, i = %d, radius[i-1] = %f \n", radius1, i, radius[(i/P_size)-1] );
        double a, b, c, d, e1, f, g;
        double deltar2 = 1.0/(grid_spacing[0]*grid_spacing[0]);
        double deltar_r = 1.0/(radius1*grid_spacing[0]);     
        double deltaphi = 1.0/(radius1*radius1*grid_spacing[1]*grid_spacing[1]);

        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r + deltaphi;
        e1 = deltaphi;
        
        
        f = (1.0/(grid_spacing[0]));
        g = -1.0*f;
        if(i<(P_size)){
            cs_entry(triplet, i, i, 1.0);
            
        }else if(i>=(P_size) && i<(2*P_size)){
            
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else if(i>=((P_size)*(P_size+1))){
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else{

            if(i%P_size==0){
                cs_entry(triplet, i, i+(P_size-1), c);
                cs_entry(triplet, i, i+(P_size-2), e1);
            }else if(i%P_size==1){
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i+(P_size-2), e1);
            }else{
                cs_entry(triplet, i, i-1, c);
                cs_entry(triplet, i, i-2, e1);
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
    cs_lusol(0, matrix, source_old, 1e-15);
    cs_spfree(triplet);
    cs_spfree(matrix);
    
    
    
    
}*/

/*
void sparse3(){
    
    const int M = (P_size*(P_size+2));
    const int N = (P_size*(P_size+2));
    int i;
    
    struct cs_sparse *triplet = cs_spalloc(M, N, 100*N, 1, 1);
    
    
    // Fill the diagonal, and the band above and below the diagonal with some
    // values.
    
    
    for (i=0; i<M; i++) {
        double radius1 = radius[i/(P_size)];
        //printf("radius = %f, i = %d, radius[i-1] = %f \n", radius1, i, radius[(i/P_size)-1] );
        double a, b, c, d, e1, f, g;
        double deltar2 = 1.0/(grid_spacing[0]*grid_spacing[0]);
        double deltar_r = 1.0/(radius1*grid_spacing[0]);     
        double deltaphi = 1.0/(radius1*radius1*grid_spacing[1]*grid_spacing[1]);

        a = deltar2;
        b = (-2.0*deltar2) - (deltar_r);
        c = -2.0*deltaphi;
        d = deltar2 + deltar_r+c;// + deltaphi;
        e1 = deltaphi;
        
        
        f = (1.0/(grid_spacing[0]));
        g = -1.0*f;
        if(i<(P_size)){
            cs_entry(triplet, i, i, 1.0);
            
        }else if(i>=(P_size) && i<(2*P_size)){
            
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else if(i>=((P_size)*(P_size+1))){
            cs_entry(triplet, i, (i-P_size), g);
            cs_entry(triplet, i, i, f);
            //printf(" d = %f, b = %f position i-P_size = %d , i = %d \n",d, b, i-P_size, i);
        }else{
            
            if(i%P_size==(P_size-1)){
                cs_entry(triplet, i, i-1,e1);// c+e1);
                cs_entry(triplet, i, i-(P_size-1), e1);
            }else if(i%P_size==0){
                cs_entry(triplet, i, i+1,e1);// c+e1);
                cs_entry(triplet, i, i+(P_size-1), e1);
            }else{
                cs_entry(triplet, i, i-1, e1);
                cs_entry(triplet, i, i+1, e1);
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
    cs_lusol(0, matrix, source_center, 1e-15);
    cs_spfree(triplet);
    cs_spfree(matrix);
    
    
    
    
}
*/


void check_pressure(){
    int i, j, position;
    boolean periodic=TRUE;
    int pertubations = (int)(P_size/pertubs);
    printf("pertubations = %d \n", pertubations);
    
    for(i=0; i<(P_size+2); i++){
        for(j=0; j<P_size; j++){
            position = (i*P_size) + j;
            
            /*if(i>0){
                if(j>P_size/2){
                    int new_pos = position-(P_size/2);
                    source[position] = source[(i*P_size) + new_pos+1];/*
                    if(source[position] != source[(i*P_size) + new_pos]){
                        if(source[position-P_size +1 + pertubations]
                    }
                } 
            }
            */
            /*if(j>=(P_size/2)){
                int new = j-(P_size/2);
                source[position] = source_old[((i+1)*P_size)-7 - new];
                //source[position- j + new] = source[position];
            }*/
            /*
            if(j>P_size/2){
                source[position] = source[position - j+(P_size/2)];
            }*/
            
            
            /*
            if(j>(P_size/2)){
                if(j>=(P_size-1-(2*pertubations))){
                    source[position] = source[(i*P_size) + (j%pertubations)];
                }/*else{
                    source[position] = source[position-(P_size/2)];

                }
            }*/
            
            /*if(j<=P_size/2){
                source[position] = source[position];
            }else{
                int new = j-(P_size/2);
                source[position] = source_old[((i+1)*P_size)-(pertubations/2)-new];

            }*/
            
            if(j>2*pertubations){
                source[position] = source[position-pertubations];
            }
            
            
            /*
            
            if(source_old[position-4] != source[position]){
                if(j<P_size/2){
                    source[position] = source[position];
                }else if(j>P_size/2){
                    source[position] = source_old[position-4];
                }else{
                    //source[position] = source_center[position-4];
                }
            }*/

        }
    }
    
}


void poisson_pressure(double *d1uphi2, double *d1ur2, double *pressure2, int P_size2, double r12, double r22, double Wi1, double Wo1, double Re2, double n2, double pertubs2)
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
    Re = Re2;
    Wi = Wi1;
    Wo = Wo1;
    r1 = r12;
    r2 = r22;
    pertubs = pertubs2;
    int i;
    double *d1uphi = &d1uphi2[0];
    double *d1ur = &d1ur2[0];
    double *pressure = &pressure2[0];

    /*Wi = 5.0;
    Wo = 10.0;
    r1 = 1.0;
    r2 = 2.0;
    

    
    d1uphi = (double*) malloc(P_size*P_size*sizeof(double));
    
    d1ur = (double*) malloc(P_size*P_size*sizeof(double));*/


    source = (double*) malloc((P_size)*(P_size+2)*sizeof(double));
    //source_old = (double*) malloc((P_size)*(P_size+2)*sizeof(double));
    //source_center = (double*) malloc((P_size)*(P_size+2)*sizeof(double));
    
    Pinner = 1.0; 
    radius = (double*) malloc((P_size+2)*sizeof(double)); 
    uphi_new = (double*) malloc(P_size*sizeof(double));

    //open_file();
    
    
    
    grid_spacing[0] = (r2-r1)/n2;
    grid_spacing[1] = (2.0*Pi)/n2;
    //printf("dr = %f, dphi = %f \n", grid_spacing[0], grid_spacing[1]);
    set_radius();
    
    fill_source(d1uphi, d1ur);
    sparse();
    //sparse2();
    //sparse3();
    check_pressure();

    
    
    // Print the solution vector.
    /*output=fopen("pressure.txt", "w");
    output1 = fopen("source.txt", "w");
    
    for (i=0; i<((P_size)*(P_size+2)); i++) {
        fprintf(output, "%16.12e ",source[i]);
        fprintf(output1, "%16.12e ",source_old[i]);
        //pressure[i] = source[i];
        
        if((i+1)%(P_size)==0){
            fprintf(output, "\n");
            fprintf(output1, "\n");
            //fprintf(output, "%f \n", source[i]);
        }
        
    }*/
    
    for (i=0; i<((P_size)*(P_size+2)); i++) {
        //printf("source[%d] = %+5.4e\n", i, source[i]);
        pressure[i] = source[i];
        
    }
    
    /*
    finite_difference();
    output3=fopen("uphi_new.txt", "w");
    
    for (i=0; i<(P_size); i++) {
        
        fprintf(output3, "%f \n", uphi_new[i]);
    }
    
    fclose(output3);
    fclose(output1);
    
    fclose(output);*/

    // Clean up memory usage.
    free(uphi_new);
    free(radius);
    //free(source_center);
    //free(source_old);
    free(source);
    
    //return 0;
}


