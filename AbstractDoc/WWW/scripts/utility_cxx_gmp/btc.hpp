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

mpf_class a = 0.9999999999999999999997_mpf;
mpf_class b = 0.9999999999999999996_mpf;

#define NUM_POINTS 1000
#define NUM_COMMANDS 3
#define PRECISION 1000

// Conversion to double: mpf_t / 10^CUT_OFF
#define CUT_OFF 12

mpf_class Catalan(mpf_t);
mpf_class def(mpf_class);
mpf_class psi(mpf_class);
mpf_class phi(mpf_class);
mpf_class always_fork(mpf_class);
unsigned long int cut_to_long_int(mpf_class);

void generate_plot(double[NUM_POINTS], double[NUM_POINTS], double[NUM_POINTS]);



