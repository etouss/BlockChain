#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
#include <fstream>
#include <string>
#include "parameters.hpp"

mpf_class a = 0.999_mpf;
mpf_class b = 0.9996_mpf;

// Number of points to plot
#define NUM_POINTS 1000

// Don't touch this. It serves to set up the pipe c++ -> terminal gnu
#define NUM_COMMANDS 3

// Precision of floats, in bits. 
#define PRECISION 1000

// Conversion to double: mpf_t / 10^CUT_OFF
// Tune this if graphs look funny. 
// This is not our fault. A normalisation factor (1-a*b) should solve this problem.
//
// Ladderish --> Diminish CUT_OFF (negative values allowed)
// Chaotic/Periodic --> Augment CUT_OFF
#define CUT_OFF 0



// Set of functions (defined in btc.cpp)
mpf_class Catalan(mpf_t);
mpf_class def(mpf_class);
mpf_class psi(mpf_class);
mpf_class phi(mpf_class);
mpf_class always_fork(mpf_class);
unsigned long int cut_to_long_int(mpf_class);

void generate_plot(double[NUM_POINTS], double[NUM_POINTS], double[NUM_POINTS]);



