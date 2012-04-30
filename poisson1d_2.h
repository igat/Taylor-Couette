/* making a header file for poisson1d_2.c */
#ifndef POISSON1D_2_H_INCLUDED
#define POISSON1D_2_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cs.h"

FILE *output;
FILE *output1;
FILE *output2;
FILE *output3;

#define Pi 3.14159265

double deriv_r(int position);

double deriv_phi(int position);

double delta_r_ur(int position);

void set_radius();

void fill_source();

void finite_difference();

void sparse();

void poisson_pressure(double *d1uphi, double *d1ur, double *pressure, int P_size, double r1, double r2, double Wi, double Wo);

#endif
