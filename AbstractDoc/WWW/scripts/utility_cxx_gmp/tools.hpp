#ifndef __TOOLS__HPP
#define __TOOLS__HPP 

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
#include <fstream>
#include <string>

// Alpha and Beta
extern mpf_class a; 
extern mpf_class b; 

// Number of points to plot
#define NUM_POINTS 500

// Don't touch this. It serves to set up the pipe c++ -> terminal gnuplot
#define NUM_COMMANDS 6

// Precision of floats, in bits. 
#define PREC 1000

// Conversion to double: mpf_t / 10^CUT_OFF
// Tune this if graphs look funny. 
// This is not our fault. A normalisation factor (1-a*b) should solve this problem.
//
// Ladderish --> Diminish CUT_OFF (negative values allowed)
// Chaotic/Periodic --> Augment CUT_OFF
#define CUT_OFF -15

// In case one has an infinite series without (implemented) closed form
#define SUM_LIMIT 150

// Set of functions (defined in tools.cpp)
mpf_class Cat(mpf_t);

// Default 
mpf_class def(mpf_class);

// Always fork
mpf_class psi(mpf_class);
mpf_class phi(mpf_class);
mpf_class always_fork(mpf_class);
// Window fork
mpf_class phi_j(int k, int j, mpf_class h);
mpf_class psi_j(int k, int j, mpf_class h);
mpf_class PHI_k(int k, mpf_class h);
mpf_class PSI_k(int k, mpf_class h);
mpf_class PHI_k_closed(int k, mpf_class h);
mpf_class PSI_k_closed(int k, mpf_class h);

mpf_class window_fork(int k, mpf_class);


unsigned long int cut_to_long_int(mpf_class);

void generate_plot(int, double[NUM_POINTS], double [][NUM_POINTS], std::string * );

// base^exp
mpf_class custom_pow(mpf_class, int);



#endif