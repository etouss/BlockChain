#include "btc.hpp"

using namespace std;


int main(int argc, char** argv){

    // Floating point precision
    mpf_set_default_prec (PRECISION);

    double h_step = (double)1/NUM_POINTS;
    mpf_class x, y;
    double h;

    // Generate each point
    double xvals[NUM_POINTS], yvals[NUM_POINTS], yvals2[NUM_POINTS];
    for(int i = 1; i < NUM_POINTS; i++){
        
        h = i*h_step;
        x = def(h);
        y = always_fork(h);
        
        xvals[i] = h;
        yvals[i] = cut_to_long_int(x);
        yvals2[i] = cut_to_long_int(y);

        gmp_printf("%f,", h);
        gmp_printf("%.*Ff\n", PRECISION/3, x);
        cout << endl << h << "," << cut_to_long_int(x) << endl;
        
    }
    generate_plot(xvals, yvals, yvals2);

    return 0;
}

// Catalan numbers generating function: (1-sqrt(1-4x))/(2x)
mpf_class Catalan(mpf_class x){

    mpf_class res;
    res = 1 - 4 * x;

    if(cmp(res,0)<0){
        cout << "Negative sqrt" << endl;
        cout << "Fatal error. Exiting. " << endl;
        exit(0);
    }

    res = (1 - sqrt(res)) / (2 * x);

    return res;
}

// Default strategy: abh/(1-a) 
mpf_class def(mpf_class h){

    mpf_class num, den, res;

    num = a * b * h;
    den = (1 - a * b);

    res = num / den;

    return res;
}


mpf_class always_fork(mpf_class h){
    mpf_class res;
    res = (1 - b) * phi(h) / (1 - psi(h));

    return res;
}


// Auxiliary function for always_fork()
mpf_class psi(mpf_class h){
    mpf_class res, fact, arg, cat;

    
    arg = (1 - h) * h * b * b * a;
    cat = Catalan(arg);
    fact = a * b * h;
    res = fact * cat;

    return res;
}

// Auxiliary function for always_fork()
mpf_class phi(mpf_class h){

    mpf_class res, b2h1mh, fact, num, den, arg, cat1, cat2;
    
    // Factor:
    num = a * b * h;
    den = (1 - a) * (1 - b);
    fact = num / den;

    // First Catalan:
    b2h1mh = b * b * h *(1 - h);
    cat1 = Catalan(b2h1mh);

    // Second Catalan * a:
    arg = a * b2h1mh;
    cat2 = Catalan(arg);

    // Result:
    res = (cat1 - a * cat2) * fact ;

    return res;
}


void generate_plot(double xvals[NUM_POINTS], double yvals[NUM_POINTS],double yvals2[NUM_POINTS]){

    string commandsForGnuplot[] = {
        "set title \"Utilities\"", 
        "plot 'Default' with lines, 'Always Fork' with lines",
        "set xtics 10"
    };
    FILE * temp1 = fopen("Default", "w");
    FILE * temp2 = fopen("Always Fork", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

    int i;
    for (i=0; i < NUM_POINTS; i++)
    {
        fprintf(temp1, "%lf %lf \n", xvals[i], yvals[i]);
        fprintf(temp2, "%lf %lf \n", xvals[i], yvals2[i]);
    }

    for (i=0; i < NUM_COMMANDS; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i].c_str()); 
    }

}

/*******************
Tune CUT_OFF in btc.hpp.
This function takes a mpf_class and translate it into a proper long int. 

If CUT_OFF == 0 :

    87521.87621387612 -- > 87522

IF CUT_OFF = 1 
 
    87521.87621387612 -- > 8752

IF CUT_OFF = -1
    
    87521.87621387612 -- > 875219

*******************/
unsigned long int cut_to_long_int(mpf_class x){

    if(!CUT_OFF) return x.get_si();

    mpf_t power_of_ten, res, mpf_t_x;
    mpf_inits(power_of_ten, res, mpf_t_x,NULL);
    mpz_t power_of_ten_z;
    mpz_init(power_of_ten_z);
    mpf_set(mpf_t_x, x.get_mpf_t());
    
    mpz_ui_pow_ui(power_of_ten_z, 10, abs(CUT_OFF));
    mpf_set_z(power_of_ten, power_of_ten_z);
    if(CUT_OFF > 0) mpf_div(res, mpf_t_x, power_of_ten);
    else mpf_mul(res, mpf_t_x, power_of_ten);
    return mpf_get_si (res);
}