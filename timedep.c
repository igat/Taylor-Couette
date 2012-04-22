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



#define Pi (4.0*atan(1.0))
#define max_iterations 50
#define max_iterations1 10


FILE *output;
FILE *output2;
FILE *output3;

static double r1, r2;
static int N[2];
static double phi_1, phi_2;
static double *U_R, *U_PHI;
static double *d1uphi;
//static double *Lr, *Lphi;
static double grid_spacing[2];
static double *radius;
static double Re;
static double V_phi_outer, V_phi_inner;
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
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            /*if(i==0){
                U_PHI[position] = V_phi_inner;
            }else if(i==(N[0]-1)){
                U_PHI[position] = V_phi_outer;
            }else{
                U_PHI[position] = 0.0;
            }*/
            U_PHI[position] = d1uphi[i];
            U_R[position] = 0.0;
            //printf("Psi[%d] = %f\n", position, Psi0[position]);
            //printf("omega[%d] = %f\n", position, Omega0[position]);
            //printf("uphi[i-1] = %f, uphi[i+1] = %f \n",u_phi[i-1], u_phi[i+1]);

        }
    
    }
    printf("Initialized psi and omega \n");
                 
            
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
        value = (-(2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]);
    }else if(position>=((N[0]-1)*N[1])){        //reflecting boundary
        //printf("at boundary condition (greater than n0-1*n1) for r2. position = %d \n", position);
        value = (-(2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]);

    }else{
        value = (w[position + N[1]] + w[position-N[1]] - (2.0*w[position]))/(grid_spacing[0]*grid_spacing[0]); 
    }
    return value;
}

double delta_phi(double *w, int position, int r_or_phi){
    double value;
    if(position%N[1]==0){
        value = (w[position + 1] - w[position + N[1]-1])/(grid_spacing[1]);
    }else if(position%N[1] ==(N[1]-1)){
        value = (w[position-N[1]+1] - w[position-1])/(grid_spacing[1]);
    }else{
        value = (w[position + 1] - w[position-1])/(grid_spacing[1]);
    }
    return value;
}

double delta_r(double *w, int position, int r_or_phi){
    double value;
    if(position<=(N[1]-1)){     //reflecting boundary
        value = (2.0*w[position + N[1]])/(grid_spacing[0]);

        //value = (2.0*w[position + N[1]])/(grid_spacing[0]);
    }else if(position>=((N[0]-1)*N[1])){        //reflecting boundary everywhere
        value = ( - 2.0*w[position-N[1]])/(grid_spacing[0]); 

        //value = (2.0*w[position-N[1]])/(grid_spacing[0]); 
        //value = 0.0;
    }else{
        value = (w[position + N[1]] - w[position-N[1]])/(grid_spacing[0]); 
    }
    return value;
}

double laplace(double *w, int position, int r_or_phi, int r_pos){
    double value;
    value = (delta_r(w, position, r_or_phi)/radius[r_pos])+ delta_r2(w, position, r_or_phi) + (delta_phi2(w, position, r_or_phi)/(radius[r_pos]*radius[r_pos]));
    return value;
}
void save_data(){
    output=fopen("u_phi32.txt", "w");
    int i, j, position;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(j==(N[1]/2)){
                fprintf(output, "%f  %f  \n", radius[i], U_PHI[position]);
            }
        }
        
    }
    fclose(output);
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
        
        printf("%f, i=%d \n", d1uphi[i], i);
        i+=1;
    }
    fclose(output2);
    
}


int main(int argc, char **argv)
{
    CFL = 0.025;
    r1 = 1.0;
    r2 = 2.0;
    N[0] = 32; // array size in each direction, N[0] = rdim
    N[1] = 32; //N[1] = PhiDim
    phi_1 = 0.0;
    phi_2 = 2.0*Pi;
    //phis go from phi = [0, 2pi]
    grid_spacing[0] =(r2-r1)/N[0];
    grid_spacing[1] = (phi_2 - phi_1)/N[1]; // grid spacing
    U_R  = (double*) malloc(N[0]*N[1]*sizeof(double));
    U_PHI  = (double*) malloc(N[0]*N[1]*sizeof(double));
    d1uphi = (double*) malloc(N[0]*sizeof(double));
    //double turn_omega = 2.0;
    open_file();

    
    V_phi_inner = 10.0;
    V_phi_outer = 0.0;
    time = 0.0;
    /* --fix these!!!
     Boundary conditions:
          
     */
    
    radius = (double*) malloc((N[0]+1)*sizeof(double));
    
    initialize_u();
    //visual_set_texdata(Psi0);
    //now we have to use a finite differencing scheme and integrate
    //the vorticity in time using RK3 and then use a relaxing scheme
    //to get psi at each step.
    

    
    const double nu = 0.1; //viscocity
    //Re = V_phi_outer*grid_spacing[1]/nu; //Reynolds number
    Re = 1.0;
    double rho;
    rho = cos(Pi/N[0]);
    omega= 2.0/(1.0+sqrt(1.0-(rho*rho)));
    printf("Re = %f, omega = %f \n", Re, omega);
    resid = 0;

    visual_init(argc, argv);
    visual_launch();


    free(radius);

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
    TextureData = (GLfloat*) malloc(N[0]*N[1]*3*sizeof(GLfloat));
    counter = 0;
}

void visual_set_texdata(double *z)
{
    double maxval = -1e13;
    double minval = +1e13;
    int n, i, j;
    for (n=0; n<(N[0]*N[1]); ++n) {
        maxval = (z[n] > maxval) ? z[n] : maxval;
        minval = (z[n] < minval) ? z[n] : minval;
    }
    for (i=0; i<N[0]; i++) {
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
    printf("texture data max/min = (%e, %e)\n", maxval, minval);
    
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
    }
    else if (key == 'p') {
        TextureMode = 'p';
    }
    else if (key == 's') {
        save_data();
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
    double L1r, L2r, L1phi, L2phi;
    
    tSTEP = CFL*grid_spacing[0]/V_phi_inner;
    //V_phi_inner += acceleration*tSTEP;
    int i, j, position;
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(i==0){
                ur1[position] = 0.0;
                up1[position] = V_phi_inner;
            }else if(i==(N[0]-1)){
                ur1[position] = 0.0;
                up1[position] = V_phi_outer;
            }else{
                L1r = -(Re*U_PHI[position]*U_PHI[position]/radius[i])+laplace(U_R, position,1, i) - (U_R[position]/(radius[i]*radius[i])) - (2.0*delta_phi(U_PHI, position, 2)/(radius[i]*radius[i]));
                L2r = Re*(-(U_PHI[position]*U_PHI[position]/radius[i]) + (U_R[position]*delta_r(U_R, position, 1)) + (U_PHI[position]*delta_phi(U_R, position, 1)/radius[i]));
                ur1[position] = U_R[position] + tSTEP*(L1r - L2r);
                
                L1phi = laplace(U_PHI, position, 2, i) + (2.0*delta_phi(U_R, position, 1)/(radius[i]*radius[i])) - (U_PHI[position]/(radius[i]*radius[i]));
                L2phi = Re*((U_R[position]*delta_r(U_PHI, position, 2)) + (U_PHI[position]*delta_phi(U_PHI, position, 2)/radius[i]) + (U_PHI[position]*U_R[position]/radius[i]));
                
                
                up1[position] = U_PHI[position] + tSTEP*(L1phi - L2phi);
            }
            
            
            //printf("L1 = %f, L2 = %f \n",L1, L2);
        }
    }
    //second integration now
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(i==0){
                ur2[position] = 0.0;
                up2[position] = V_phi_inner;
            }else if(i==(N[0]-1)){
                ur2[position] = 0.0;
                up2[position] = V_phi_outer;
            }else{
                L1r = -(Re*up1[position]*up1[position]/radius[i])+laplace(ur1, position,1, i) - (ur1[position]/(radius[i]*radius[i])) - (2.0*delta_phi(up1, position, 2)/(radius[i]*radius[i]));
                L2r = Re*(-(up1[position]*up1[position]/radius[i]) + (ur1[position]*delta_r(ur1, position, 1)) + (up1[position]*delta_phi(ur1, position, 1)/radius[i]));
                ur2[position] = ((3.0/4.0)*U_R[position]) + (ur1[position]/4.0) + (tSTEP*(L1r - L2r)/4.0);
                
                L1phi =laplace(up1, position, 2, i) + (2.0*delta_phi(ur1, position, 1)/(radius[i]*radius[i])) - (up1[position]/(radius[i]*radius[i]));
                L2phi = Re*((ur1[position]*delta_r(up1, position, 2)) + (up1[position]*delta_phi(up1, position, 2)/radius[i]) + (up1[position]*ur1[position]/radius[i]));
                
                up2[position] = ((3.0/4.0)*U_PHI[position]) + (up1[position]/4.0) + (tSTEP*(L1phi - L2phi)/4.0);
            }
            
            
            
        }
    }
    //third integration
    for(i=0; i<N[0]; i++){
        for(j=0; j<N[1]; j++){
            position = (i*N[1]) + j;
            if(i==0){
                ur3[position] = 0.0;
                up3[position] = V_phi_inner;
            }else if(i==(N[0]-1)){
                ur3[position] = 0.0;
                up3[position] = V_phi_outer;
            }else{
                L1r = -(Re*up2[position]*up2[position]/radius[i])+laplace(ur2, position,1, i) - (ur2[position]/(radius[i]*radius[i])) - (2.0*delta_phi(up2, position, 2)/(radius[i]*radius[i]));
                L2r = Re*(-(up2[position]*up2[position]/radius[i]) + (ur2[position]*delta_r(ur2, position, 1)) + (up2[position]*delta_phi(ur2, position, 1)/radius[i]));
                ur3[position] = ((1.0/3.0)*U_R[position]) + (2.0*ur2[position]/3.0) + (2.0*tSTEP*(L1r - L2r)/3.0);
                
                L1phi = laplace(up2, position, 2, i) + (2.0*delta_phi(ur2, position, 1)/(radius[i]*radius[i])) - (up2[position]/(radius[i]*radius[i]));
                L2phi = Re*((ur2[position]*delta_r(up2, position, 2)) + (up2[position]*delta_phi(up2, position, 2)/radius[i]) + (up2[position]*ur2[position]/radius[i]));
                
                up3[position] = ((1.0/3.0)*U_PHI[position]) + (2.0*up2[position]/3.0) + (2.0*tSTEP*(L1phi - L2phi)/3.0);
            }
            
            
            
        }
    }
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
    if(difference<=(1e-8)){
        save_data();
        glutDestroyWindow(WindowID);
        exit(0);
        
    }
    
}
