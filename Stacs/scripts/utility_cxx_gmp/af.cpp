#include "af.hpp"

using namespace std;

int main(int argc, char** argv){

    double left_h = 0;
    double right_h = 1;

    if(argc == 3){
        double arg1 = atof(argv[1]);
        double arg2 = atof(argv[2]);
        if((arg1 >= 0 && arg1 <= 1) &&
           (arg2 >= 0 && arg2 <= 1) &&
           (arg2  > arg1))
        {
            left_h = arg1;
            right_h = arg2;
        }
        else{
            cout << "Please, be serious. " << endl;
            exit(0);
        }
    }
    
    generate_af_plot(left_h, right_h);
    return 0;
}


void generate_af_plot(double initial_i, double final_i){

    mpf_set_default_prec (PREC);
    double h_step = (double)(final_i - initial_i)/NUM_POINTS;
    double h;


    // Generate each point
    double xvals[NUM_POINTS];

    double Y[2][NUM_POINTS];
    
    for(int i = 0; i < NUM_POINTS; i++){
        h = initial_i + (i+1) * h_step;
        // Gnuplots' fault:
        if(i == NUM_POINTS-1){ h -= h_step;}

        xvals[i] = h;
        Y[0][i] = cut_to_long_int(def(h));
        Y[1][i] = cut_to_long_int(always_fork(h));
        

        /*gmp_printf("%f,", h);
        gmp_printf("%.*Ff\n", PREC/3, Y[0][i]);
        cout << endl << h << "," << cut_to_long_int(def(h)) << endl;*/
        cout << i << endl;
    }

    string names[2] = {"Default.temp", "Always_Fork.temp"};
    
    generate_plot(2, xvals, Y, names);
}