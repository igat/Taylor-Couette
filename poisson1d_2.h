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

double deriv_r(int position, double *d1uphi);

double deriv_phi(int position, double *d1ur1);

double delta_r_ur(int position, double *d1ur1);

void set_radius();

void deriv_phi_new(double *new_uphi, double *old_uphi);

double deriv_phi2(double *w, int position);


void fill_source(double *d1uphi, double *d1ur);

void finite_difference();

void sparse();

void poisson_pressure(double *d1uphi2, double *d1ur2, double *pressure, int P_size, double r1, double r2, double Wi, double Wo, double Re2, double n2, double pertubs);

#endif
