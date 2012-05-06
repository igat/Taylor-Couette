/* to run:
 
 gcc -Wall -o timedep timedep.c -L/usr/X11/lib -lglut -lGLU -lGL -lXmu -lXext -lXi -lX11 -I/usr/X11/include

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <GL/glut.h>


/*
 
 This is a code that uses the 2D cylindrical coordinates of the 
 Navier-Stokes Equations. It assumes there is zero gradient
 in the Pressure and body forces (or assume they are zero)
 and that the fluid is incompressible.
 
*/



#define Pi 3.1415926535897932384626433832795028841971693993751058209
#define max_iterations 500
#define max_iterations1 10


FILE *output;
FILE *output2;
FILE *output3;

static double r1, r2;
static int N[2];
static double phi_1, phi_2;
static double *U_R, *U_PHI, *Pressure;
static double *d1uphi;
//static double *Lr, *Lphi;
static double grid_spacing[2];
static double *radius, *phi_value;
static double Re;
static double V_phi_outer, V_phi_inner, Wi, Wo;
static double omega;
static double residual[max_iterations];
static int resid;
static double CFL;
static double time, tSTEP;
static char TextureMode = 'p';

void visual_init(int argc, char **argv);
void visual_set_texdata(double *z);
void visual_launch();

void initialize_u(){
    int i, j, position;
    for(i=0; i<N[0]; i++){
        radius[i] = r1 + (i*grid_spacing[0]);
        phi_value[i] = phi_1 + (i*grid_spacing[1]);

        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            U_PHI[position] = d1uphi[i];
            U_R[position] = 0.0;

        }
    
    }
    //not quite sure if this next part is right
    for(i=0; i<(N[0]+2); i++){
        //radius[i] = r1 + (i*grid_spacing[0]);
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            Pressure[position] = 1.0;
        }
        
    }
    
    printf("Initialized pressure and velocities \n");
                 
            
}

double delta_phi2(double *w, int position, int r_or_phi){
    double value;
    if(position%N[1]==0){
        //printf("at boundary condition (0) for phi2. position = %d \n", position);
        value = (w[position + 1] + w[position + N[1]-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else if(position%N[1] ==(N[1]-1)){
        //printf("at boundary condition (n1-1) for phi2. position = %d \n", position);
        value = (w[position-N[1]+1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }else{
        value = (w[position + 1] + w[position-1] - (2.0*w[position]))/(grid_spacing[1]*grid_spacing[1]);
    }
    return value;
}

double delta_r2(double *w, int position, int r_or_phi){
    double value;
    if(position<=(N[1]-1)){     //reflecting boundary
        if(r_or_phi==1){
            value = (-(2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]);
        }else{
            value = 0.0;
        }
    }else if(position>=((N[0]-1)*N[1])){        //reflecting boundary
        //printf("at boundary condition (greater than n0-1*n1) for r2. position = %d \n", position);
        if(r_or_phi==1){
            value = (w[position-N[1]]-(2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]);
        }else{
            value = 0.0;
        }
        
    }else{
        value = (w[position + N[1]] + w[position-N[1]] - (2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]); 
        //value = (w[position] + w[position-(2*N[1])] - (2.0*w[position-N[1]]))/(grid_spacing[0]*grid_spacing[0]); 
    }
    return value;
}

double delta_phi(double *w, int position, int r_or_phi){
    double value;
    if(position%N[1]==0){
        value = (w[position + 1] - w[position + N[1]-1])/(2.0*grid_spacing[1]);
    }else if(position%N[1] ==(N[1]-1)){
        value = (w[position-N[1]+1] - w[position-1])/(2.0*grid_spacing[1]);
    }else{
        value = (w[position + 1] - w[position-1])/(2.0*grid_spacing[1]);
    }
    return value;
}

double delta_r(double *w, int position, int r_or_phi){
    double value;
    if(position<=(N[1]-1)){     //reflecting boundary
        if(r_or_phi==1){
            value = (2.0*w[position + N[1]])/(2.0*grid_spacing[0]);
        }else{
            value = Wi;
        }
        
        //value = (2.0*w[position + N[1]])/(grid_spacing[0]);
    }else if(position>=((N[0]-1)*N[1])){        //reflecting boundary everywhere
        if(r_or_phi==1){
            value = (-w[position-N[1]])/(2.0*grid_spacing[0]); 
        }else{
            value = Wo;
        }
        
        //value = (2.0*w[position-N[1]])/(grid_spacing[0]); 
        //value = 0.0;
    }else{
        //value = (w[position + N[1]] - w[position-N[1]])/(grid_spacing[0]); 
        value = (w[position+N[1]] - w[position-N[1]])/(2.0*grid_spacing[0]); 
    }
    return value;
}

double laplace(double *w, int position, int r_or_phi, int r_pos){
    double value;
    value = (delta_r(w, position, r_or_phi)/radius[r_pos])+ delta_r2(w, position, r_or_phi) + (delta_phi2(w, position, r_or_phi)/(radius[r_pos]*radius[r_pos]));
    return value;
}

void save_data(double *ur, double *uphi){
    output=fopen("Pressure_timedep2.txt", "w");
    output2=fopen("ur_timedep2.txt", "w");
    output3 = fopen("uphi_timedep2.txt", "w");
    int i, j, position;
    for(i=0; i<(N[0]+2); i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            fprintf(output, "%16.12e  ", Pressure[position]);
            if(j==(N[1]-1)){
                fprintf(output, "\n");
                //fprintf(output, "%f  %f  \n", radius[i], U_PHI[position]);
            }
        }
        
    }
    
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            fprintf(output2, "%16.12e  ", ur[position]);
            fprintf(output3, "%16.12e  ", uphi[position]);
            
            if(j==(N[1]-1)){
                fprintf(output2, "\n");
                fprintf(output3, "\n");
                //fprintf(output, "%f  %f  \n", radius[i], U_PHI[position]);
            }
        }
        
    }
    fclose(output3);
    fclose(output2);
    fclose(output);
}


double delta_Pressure_r(int position){
    double value;
    value = (Pressure[position+N[0]]-Pressure[position])/(grid_spacing[0]);
    //printf("Pressure derivative = %f, uphi[%d]^2/radius = %f \n",value,position, U_PHI[position]*U_PHI[position]/radius[position/N[0]]);
    return value;
}

double delta_Pressure_phi(int position){
    double value;
    if(position%N[1]==0){
        value = (Pressure[position+N[1]+1] - Pressure[position + (2*N[1]) - 1])/(2.0*grid_spacing[1]);
    }else if(position%N[1] ==(N[1]-1)){
        value = (Pressure[position+1] - Pressure[position + N[1] - 1])/(2.0*grid_spacing[1]);
    }else{
        value = (Pressure[position+N[1]+1] - Pressure[position + N[1] - 1])/(2.0*grid_spacing[1]);
    }
    //printf("Pressure derivative = %E, uphi[%d]^2/radius = %f \n",value,position, U_PHI[position]*U_PHI[position]/radius[position/N[0]]);
    
    return value;
}

double deriv_r(int position, double *d1uphi){
    double value;
    value = (d1uphi[position] - d1uphi[position-N[0]])/(grid_spacing[0]); 
    //printf("derivative in r of uphi = %f, d1uphi[%d] = %f \n", value, position, d1uphi[position]);
    return value;
}

double deriv_phi(int position, double *d1ur1){
    double value;
    double *d1ur  = &d1ur1[0];
    if(position%N[0]==0){
        value = (d1ur[position+1] - d1ur[position+N[0] - 1])/(2.0*grid_spacing[1]);
    }if((position+1)%N[0]==0){
        //printf("here, position = %d, new position = %d \n", position, position-P_size +1);
        value = (d1ur[position-N[0] + 1] - d1ur[position - 1])/(2.0*grid_spacing[1]);
    }else{
        value = (d1ur[position+1] - d1ur[position-1])/(2.0*grid_spacing[1]); 
    }
    
    //printf("derivative in phi of ur = %f, d1ur[%d] = %f \n", value, position, d1ur[position]);
    return value;
}

double delta_r_ur(int position, double *d1ur1){
    double value;
    double *d1ur  = &d1ur1[0];
    if(position>=(N[0]-1)*N[1]){
        value = (-d1ur[position-N[0]])/(grid_spacing[0]);
    }else{
        value = (d1ur[position] - d1ur[position-N[0]])/(grid_spacing[0]);
    } 
    //printf("derivative of ur in r = %E, d1ur[%d] = %E \n", value, position, d1ur[position]);
    return value;
}

void deriv_phi_first(double *value_new, double *d1uphi){
    int i, j, position2;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position2 = (i*N[1]) + j;
            value_new[position2] = delta_phi(d1uphi, position2, 2);
        }
    }
}

double delta_r_phi(int position, double *d1uphi){
    double *value_new = (double*) malloc(N[0]*N[1]*sizeof(double));
    double value;
    
    int i, j, position2;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position2 = (i*N[1]) + j;
            value_new[position2] = delta_phi(d1uphi, position2, 2);
        }
    }
    value = deriv_r(position, value_new);

    free(value_new);
    return value;
}

void relax_pressure(double *ur, double *uphi){
    const double EPS = 2.0e-16;
    double resid_sum;
    
    double *value_new = (double*) malloc(N[0]*N[1]*sizeof(double));
    deriv_phi_first(value_new, uphi);
    
    for(resid=0; resid<max_iterations; resid++){
        //double *pressure_new = (double*) malloc(N[0]*N[1]*sizeof(double));
        double residual1, source;
        double Piminus2, Piminus1, Pjplus1, Pjminus1;
        int i, j, position;
        double P_new;
        double a, b, c, d, f, g;
        resid_sum = 0.0;
        for(i=0; i<(N[0]+2); i++){
            for(j=0; j<N[1]; j++){
                position = (i*N[1]) + j;

                if(i==0){
                    P_new = 1.0;
                }else if(i==1){
                    source = r1*Wi*Wi;
                    residual1 = Pressure[position] - Pressure[position-N[1]] - (grid_spacing[0]*source/omega);
                    P_new = Pressure[position] - (omega*residual1);
                    resid_sum+=fabs(residual1);
                }else if(i==(N[0]+1)){
                    source = r2*Wo*Wo;
                    residual1 = Pressure[position] - Pressure[position-N[1]] - (grid_spacing[0]*source/omega);
                    P_new = Pressure[position] - (omega*residual1);
                    resid_sum+=fabs(residual1);
                }else{
                    //source = (2.0)*(((deriv_r((position-N[1]), uphi)/radius[i-1])*(uphi[(position-N[1])] - deriv_phi((position-N[1]), ur))) - (delta_r_ur((position-N[1]), ur)*delta_r_ur((position-N[1]), ur)) - ((1.0/(Re*radius[i-1]*radius[i-1]))*(delta_r_ur((position-N[1]), ur) + (deriv_r((position-N[1]), value_new)) + (delta_phi2(ur, (position - N[1]), 1)/radius[i-1]))));
                    source = (2.0)*(((delta_r(uphi,(position-N[1]), 2)/radius[i-1])*(uphi[(position-N[1])] - delta_phi(ur, (position-N[1]), 1))) - (delta_r(ur, (position-N[1]), 1)*delta_r(ur,(position-N[1]), 1)) - ((1.0/(Re*radius[i-1]*radius[i-1]))*(delta_r(ur, (position-N[1]), 1) + (delta_r(value_new,(position-N[1]), 2)) + (delta_phi2(ur, (position - N[1]), 1)/radius[i-1]))));
                    
                    
                    a = grid_spacing[1]*grid_spacing[1]*radius[i-1]*radius[i-1];
                    b = (-2.0*a) - (grid_spacing[0]*grid_spacing[1]*grid_spacing[1]*radius[i-1]);
                    c = grid_spacing[0]*grid_spacing[0];
                    d = c;
                    g = a + (grid_spacing[1]*grid_spacing[1]*grid_spacing[0]*radius[i-1]) - (2.0*c);
                    f = grid_spacing[0]*grid_spacing[0]*grid_spacing[1]*grid_spacing[1]*radius[i-1]*radius[i-1];
                    if(j==0){
                        Pjminus1 = Pressure[position + N[1] -1];
                        Pjplus1 = Pressure[position + 1];
                    }else if(j==(N[1]-1)){
                        Pjminus1 = Pressure[position-1];
                        Pjplus1 = Pressure[position + 1 - N[1]];
                    }else{
                        Pjminus1 = Pressure[position-1];
                        Pjplus1 = Pressure[position + 1];
                    }
                    Piminus1 = Pressure[position - N[1]];
                    Piminus2 = Pressure[position - (2*N[1])];
                    
                    
                    residual1 = (Piminus2*a) + (Piminus1*b) + (Pjplus1*c) + (Pjminus1*d) + (Pressure[position]*g) - (f*source/omega);
                    
                    P_new = Pressure[position] - (omega*residual1/g);
                    resid_sum +=fabs(residual1);

                }
                
                Pressure[position] = P_new;
                //if(resid==499){printf("Pressure[%d] = %f \n", position, Pressure[position]);}

            }
        }
        //printf("residual here = %E, resid = %d \n", resid_sum, resid);
        residual[resid] = resid_sum/(N[0]*N[1]);
        //free(pressure_new);
        //if(resid=10) printf("residual here = %E, eps = %E \n", residual[resid], EPS);
        if(resid>=2){
            double value1 = residual[0]*EPS;
            //printf("value1 = %E \n", value1);
            if(residual[resid]<=value1){
                residual[max_iterations-1] = residual[resid];
                printf("here!! residual here = %E, resid = %d \n", residual[resid], resid);
                //resid = max_iterations;
                return;
                
            
            }
        }
        //visual_set_texdata(Pressure);

        //resid_sum = 0.0;
        
    }
    free(value_new);

    
    printf("residual here = %E, resid = %d\n", residual[resid-1], resid);
    printf("Check!, Pressure[10] - Pressure[9] = %f\n",Pressure[10] - Pressure[9]);
    int i, j, position;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(j==(N[1]/2)){
                //printf("Pressure[%d] = %f \n",position, Pressure[position]);
            }
        }
    }


}





void integrate_u();

/*
void integration(){
    const double EPS = 1.0e-8;
    for(resid=0; resid<max_iterations; resid++){
        relax_psi();
        if(resid>1){
            if(residual[resid]<EPS){
                resid = max_iterations;
            }
        }
    }
    //visual_set_texdata(Psi0);
}
*/

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

void integration(double *ur_new, double *up_new, double *ur_old, double *up_old, int one_two_or_three){
    int i, j, position;
    double L1r, L2r, L1phi, L2phi;
    
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(i==0){
                ur_new[position] = 0.0;
                up_new[position] = V_phi_inner;
            }else if(i==(N[0]-1)){
                L1r = -delta_Pressure_r(position)+  ((1.0/Re)*(laplace(ur_old, position,1, i) - (ur_old[position]/(radius[i]*radius[i])) - (2.0*delta_phi(up_old, position, 2)/(radius[i]*radius[i]))));
                L2r = (-(up_old[position]*up_old[position]/radius[i]) + (ur_old[position]*delta_r(ur_old, position, 1)) + (up_old[position]*delta_phi(ur_old, position, 1)/radius[i]));
                if(one_two_or_three==1){
                    ur_new[position] = ur_old[position] + tSTEP*(L1r - L2r);
                }else if(one_two_or_three==2){
                    ur_new[position] = ((3.0/4.0)*U_R[position]) + (ur_old[position]/4.0) + (tSTEP*(L1r - L2r)/4.0);
                }else{
                    ur_new[position] = ((1.0/3.0)*U_R[position]) + (2.0*ur_old[position]/3.0) + (2.0*tSTEP*(L1r - L2r)/3.0);
                }
                //ur_new[position]=0.0;

                up_new[position] = d1uphi[i];
            }else{
                L1r = -delta_Pressure_r(position)+  ((1.0/Re)*(laplace(ur_old, position,1, i) - (ur_old[position]/(radius[i]*radius[i])) - (2.0*delta_phi(up_old, position, 2)/(radius[i]*radius[i]))));
                L2r = (-(up_old[position]*up_old[position]/radius[i]) + (ur_old[position]*delta_r(ur_old, position, 1)) + (up_old[position]*delta_phi(ur_old, position, 1)/radius[i]));
                
                L1phi = -(delta_Pressure_phi(position)/radius[i]) + ((1.0/Re)*(laplace(up_old, position, 2, i) + (2.0*delta_phi(ur_old, position, 1)/(radius[i]*radius[i])) - (up_old[position]/(radius[i]*radius[i]))));
                L2phi = ((ur_old[position]*delta_r(up_old, position, 2)) + (up_old[position]*delta_phi(up_old, position, 2)/radius[i]) + (up_old[position]*ur_old[position]/radius[i]));
                
                if(one_two_or_three==1){
                    ur_new[position] = ur_old[position] + tSTEP*(L1r - L2r);
                    up_new[position] = up_old[position] + tSTEP*(L1phi - L2phi);

                }else if(one_two_or_three==2){
                    ur_new[position] = ((3.0/4.0)*U_R[position]) + (ur_old[position]/4.0) + (tSTEP*(L1r - L2r)/4.0);
                    up_new[position] = ((3.0/4.0)*U_PHI[position]) + (up_old[position]/4.0) + (tSTEP*(L1phi - L2phi)/4.0);
                }else{
                    ur_new[position] = ((1.0/3.0)*U_R[position]) + (2.0*ur_old[position]/3.0) + (2.0*tSTEP*(L1r - L2r)/3.0);
                    up_new[position] = ((1.0/3.0)*U_PHI[position]) + (2.0*up_old[position]/3.0) + (2.0*tSTEP*(L1phi - L2phi)/3.0);
                }
                
                
                
            }
            
            
            //printf("L1 = %f, L2 = %f \n",L1, L2);
        }
    }
    
    
}


int main(int argc, char **argv)
{
    //CFL = 0.07;
    CFL = 0.0006;
    r1 = 1.0;
    r2 = 2.0;
    N[0] = 100; // array size in each direction, N[0] = rdim
    N[1] = 100; //N[1] = PhiDim
    phi_1 = 0.0;
    phi_2 = 2.0*Pi;
    //phis go from phi = [0, 2pi]
    grid_spacing[0] =(r2-r1)/100.0;
    grid_spacing[1] = (phi_2 - phi_1)/100.0; // grid spacing
    U_R  = (double*) malloc(N[0]*N[1]*sizeof(double));
    U_PHI  = (double*) malloc(N[0]*N[1]*sizeof(double));
    Pressure = (double*) malloc((N[0]+2)*N[1]*sizeof(double));

    d1uphi = (double*) malloc(N[0]*sizeof(double));
    //double turn_omega = 2.0;
    open_file();

    
    Wi = 10.0;
    Wo = 0.0;
    V_phi_inner = r1*Wi;
    V_phi_outer = r2*Wo;
    time = 0.0;
    /* --fix these!!!
     Boundary conditions:
          
     */
    
    radius = (double*) malloc((N[0])*sizeof(double));
    phi_value = (double*) malloc((N[1])*sizeof(double));    
    //now we have to use a finite differencing scheme and integrate
    //the vorticity in time using RK3 and then use a relaxing scheme
    //to get psi at each step.
    
    
    const double nu = 0.1; //viscocity
    //Re = V_phi_outer*grid_spacing[1]/nu; //Reynolds number
    Re = 1.0;
    double rho;
    //rho = cos(Pi/(N[0]-97.9981552));
    //rho = cos(Pi/(N[0]-97.99));
    rho = cos(Pi/(2.00184));
    omega= 2.0/(1.0+sqrt(1.0-(rho*rho)));
    printf("Re = %f, omega = %f \n", Re, omega);
    resid = 0;
    initialize_u();

    //relax_pressure(U_R, U_PHI, 1);



    visual_init(argc, argv);
    visual_launch();


    free(phi_value);
    
    
    free(radius);
    free(d1uphi);
    free(Pressure);
    free(U_PHI);
    free(U_R);
    
    return 0;

}

static void DrawGLScene();
static void ResizeGLScene(int Width, int Height);
static void KeyPressed(unsigned char key, int x, int y);
static void SpecialKeyPressed(int key, int x, int y);
static int counter;
//static double time = 0.0;
static int WindowID;
static int WindowWidth    = 768;
static int WindowHeight   = 768;

static float xTranslate = 0.0;
static float yTranslate = 0.0;
static float zTranslate = 2.5;
static float RotationAngleX = 0;
static float RotationAngleY = 0;

static GLuint TextureMap;
static GLfloat *TextureData;


void visual_launch()
{
    glutMainLoop();
}

void visual_init(int argc, char **argv)
{
    glutInit(&argc, argv);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
    glutInitWindowSize(WindowWidth, WindowHeight);
    glutInitWindowPosition(80, 80);
    WindowID = glutCreateWindow("2D Euler");
    
    glutDisplayFunc       (DrawGLScene);
    glutIdleFunc          (DrawGLScene);
    glutReshapeFunc       (ResizeGLScene);
    glutKeyboardFunc      (KeyPressed);
    glutSpecialFunc       (SpecialKeyPressed);
    
    glClearColor(0.2, 0.2, 0.1, 0.0);
    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (float) WindowWidth / WindowHeight, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
    
    // Establish texture id's and data
    glGenTextures(1, &TextureMap);
    TextureData = (GLfloat*) malloc((N[0])*N[1]*3*sizeof(GLfloat));
    counter = 0;
}

void visual_set_texdata(double *z)
{
    double maxval = -1e13;
    double minval = +1e13;
    int n, i, j;
    for (n=0; n<((N[0])*N[1]); ++n) {
        maxval = (z[n] > maxval) ? z[n] : maxval;
        minval = (z[n] < minval) ? z[n] : minval;
    }
    for (i=0; i<(N[0]); i++) {
        for (j=0; j<N[1]; j++) {
            
            const int m = i*N[1] + j;
            
            float ColorWidth = 0.1;
            float A = pow(0.3 / ColorWidth, 2);
            float B = pow(0.3 / ColorWidth, 2);
            float C = pow(0.3 / ColorWidth, 2);
            float v = (z[m] - minval) / (maxval - minval);
            float r = exp(-A*pow(v-0.75,2));
            float g = exp(-B*pow(v-0.50,2));
            float b = exp(-C*pow(v-0.25,2));
            
            TextureData[3*m+0] = r;
            TextureData[3*m+1] = g;
            TextureData[3*m+2] = b;
        }
    }
    printf("texture data max/min = (%16.12e, %16.12e)\n", maxval, minval);
    
    glBindTexture(GL_TEXTURE_2D, TextureMap);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, N[0], N[1], 0, GL_RGB, GL_FLOAT, TextureData);
}


static void DrawGLScene()
{
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(-xTranslate, -yTranslate, -zTranslate);
    glRotatef(90, 0, 0, 1);
    
    
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, TextureMap);
    
    glBegin(GL_QUADS);
    glNormal3f(0, 0, 1);
    glTexCoord2f(0, 0); glVertex3f(-1,-1, 0);
    glTexCoord2f(0, 1); glVertex3f(-1, 1, 0);
    glTexCoord2f(1, 1); glVertex3f( 1, 1, 0);
    glTexCoord2f(1, 0); glVertex3f( 1,-1, 0);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    glFlush();
    glutSwapBuffers();
    
    integrate_u();
}


void ResizeGLScene(int Width, int Height)
{
    if (Height==0) Height = 1;
    glViewport(0, 0, Width, Height);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    gluPerspective(45.0f, (GLfloat)Width/(GLfloat)Height, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
}
void SpecialKeyPressed(int key, int x, int y)
{
    if      (key == GLUT_KEY_RIGHT) RotationAngleY += 3;
    else if (key == GLUT_KEY_LEFT ) RotationAngleY -= 3;
    
    if      (key == GLUT_KEY_DOWN ) RotationAngleX += 3;
    else if (key == GLUT_KEY_UP   ) RotationAngleX -= 3;
    
    glutPostRedisplay();
}

void KeyPressed(unsigned char key, int x, int y)
{
    const int ESCAPE = 27;
    
    if (key == ESCAPE || key == 'q') {
        glutDestroyWindow(WindowID);
        exit(0);
    }
    else if (key == 'a') {
        integrate_u();
    }
    else if (key == 'r') {
        TextureMode = 'r';
        visual_set_texdata(U_R);

    }
    else if (key == 'p') {
        TextureMode = 'p';
        visual_set_texdata(U_PHI);

    }
    else if (key == 'o') {
        TextureMode = 'o';
        visual_set_texdata(Pressure);

    }
    else if (key == 's') {
        save_data(U_R, U_PHI);
    }
    glutPostRedisplay();
}




void integrate_u(){
    double *ur1 = (double*) malloc(N[0]*N[1]*sizeof(double));
    double *ur2 = (double*) malloc(N[0]*N[1]*sizeof(double));
    double *ur3 = (double*) malloc(N[0]*N[1]*sizeof(double));
    double *up1 = (double*) malloc(N[0]*N[1]*sizeof(double));
    double *up2 = (double*) malloc(N[0]*N[1]*sizeof(double));
    double *up3 = (double*) malloc(N[0]*N[1]*sizeof(double));
    int i, j, position;

    relax_pressure(U_R, U_PHI);
    
    //save_data();

    
    //tSTEP = CFL*grid_spacing[0]*grid_spacing[1]/((V_phi_outer) + (V_phi_inner));
    tSTEP = CFL*(r2-r1)/((V_phi_outer) + (V_phi_inner));
    //tSTEP = CFL*grid_spacing[0]*grid_spacing[1];

    integration(ur1, up1, U_R, U_PHI, 1);
    //save_data(ur1, up1);
    //return;
    integration(ur2, up2, ur1, up1, 2);

    integration(ur3, up3, ur2, up2, 3);
    //return;

    //relax_pressure(ur3, up3, 2);

    double difference = 0.0;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            difference += fabs(U_PHI[position] - up3[position]);
            U_R[position] = ur3[position];
            U_PHI[position] = up3[position];
        }
    }
    free(up3);
    free(up2);
    free(up1);
    free(ur3);
    free(ur2);
    free(ur1);
    time +=tSTEP;
    printf("Current time = %f , difference = %E \n", time, difference);
    if (TextureMode == 'r') {
        visual_set_texdata(U_R);
    }
    else if (TextureMode == 'p') {
        visual_set_texdata(U_PHI);
    }
    else if  (TextureMode == 'o') {
        visual_set_texdata(Pressure);
    }
    /*if(difference<=(1e-8)){
        //save_data();
        glutDestroyWindow(WindowID);
        exit(0);
        
    }*/
    
}
