#ifndef POISSON1D_2_H_INCLUDED
#define POISSON1D_2_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define max_iterations 1000 

FILE *output;
FILE *output1;
FILE *output2;
FILE *output3;

double RK(double value, double value2);

double RK_with_update(double value, double value2, int position);

void initialize(double r1);

void find_s(double s);

void steady(double Vinner, double Vouter, int rDIM, double *u_phi1, double rinner, double router);

#endif
