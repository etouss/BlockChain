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
#define NUM_POINTS 50000

// Don't touch this. It serves to set up the pipe c++ -> terminal gnuplot
#define NUM_COMMANDS 6

// Precision of floats, in bits. 
#define PREC 10000

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

void generate_plot_mpf(int length, double xvals[NUM_POINTS], mpf_class yvals [][NUM_POINTS], std::string * names);

// base^exp
mpf_class custom_pow(mpf_class, int);

mpz_class fact(int x);
mpz_class choose(int x, int y);
mpz_class choose_gen(int i, int j);
mpf_class Cat_g(mpf_class x, int aa, int bb);
mpz_class f(int aa, int r, int bb);
mpz_class Pent_aux(int aa, int bb, int r);
mpz_class Pent(int aa, int bb, int r);


mpf_class a1(int k, int j, int l, mpf_class h);
mpf_class sum_a1(int k, int l, mpf_class h);

mpf_class a2(int k, int j, int l, mpf_class h);
mpf_class sum_a2(int k, int l, mpf_class h);

mpf_class b1(int k, int j, int l, mpf_class h);
mpf_class sum_b1(int k, int l, mpf_class h);

mpf_class b2(int k, int j, int l, mpf_class h);
mpf_class sum_b2(int k, int l, mpf_class h);

mpf_class E(int aa, int bb);

mpf_class c(int k, int j, int l, mpf_class h);
mpf_class sum_c(int k, int l, mpf_class h);

mpf_class PP(int l, int j, mpf_class x);
mpf_class d(int k, int j, int l, mpf_class h);
mpf_class sum_d(int k, int l, mpf_class h);

mpf_class util(int k, int l, mpf_class h);



#endif
