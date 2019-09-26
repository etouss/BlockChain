#include "window_fork.hpp"

using namespace std;

int main(int argc, char** argv){

    /*mpf_class x = 0.1_mpf;
    cout << PSI_k(1, x) << endl;
    cout << PSI_k_closed(1, x) << endl;
    exit(0);    */

    double left_h = 0;
    double right_h = 1;
    int n = 2;

    if(argc == 4){
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
        n = atof(argv[3]);
    }
    
    generate_window_fork_plot(left_h, right_h, n);
    return 0;
}


void generate_window_fork_plot(double initial_i, double final_i, int n){

    mpf_set_default_prec (PREC);
    double h_step = (double)(final_i - initial_i)/NUM_POINTS;
    double h;


    // Generate each point
    double xvals[NUM_POINTS];

    double Y[n+2][NUM_POINTS];

    cout << "Computing..." << endl;
        
    for(int i = 0; i < NUM_POINTS; i++){
        h = initial_i + (i+1) * h_step;
        // Gnuplots' fault:
        if(i == NUM_POINTS-1){ h -= h_step;}

        xvals[i] = h;
        Y[0][i] = cut_to_long_int(def(h));
        for(int k = 0 ; k < n; k++){
            Y[k+1][i] = cut_to_long_int(window_fork(k+1,h));
        }
        Y[n+1][i] = cut_to_long_int(always_fork(h));


        /*gmp_printf("%f,", h);
        gmp_printf("%.*Ff\n", PREC/3, Y[0][i]);
        cout << endl << h << "," << cut_to_long_int(def(h)) << endl;*/

    }

    string names[n+2];
    names[0] = "Default.temp";
    for(int k = 0 ; k < n; k++){
        names[k+1] = "Window_Fork_" + to_string(k+1)+".temp";
    }
    names[n+1] = "Always_Fork.temp";
    
    generate_plot(n+2, xvals, Y, names);
}