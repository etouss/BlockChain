#include "gup.hpp"

using namespace std;

int main(int argc, char** argv){
    

    mpf_class max = 0;
    mpf_class k_max = 0;
    mpf_class l_max = 0;
    int LIMIT = 100;
    

    //generate_gup_plot(0, 1);
    
    //return 0;

    mpf_class min = 1000000000000000000000000_mpf;
    double i_min = 0;
    double step = 0.00000001;
    mpf_class x;
    for (double i = 0.4; i < 0.6; i = i + step){
        x = abs(def(i) - always_fork(i)) ;
        if(x < min){
            i_min = i;
            min = x;
        }
    }

    cout << "4-5: " <<i_min << endl;
    
}


void generate_gup_plot(double initial_i, double final_i){

    mpf_set_default_prec (PREC);
    double h_step = (double)(final_i - initial_i)/NUM_POINTS;
    double h;


    // Generate each point
    double xvals[NUM_POINTS];

    mpf_class Y[5][NUM_POINTS];
    
    for(int i = 0; i < NUM_POINTS; i++){
        h = initial_i + (i+1) * h_step;
        // Gnuplots' fault:
        if(i == NUM_POINTS-1){ h -= h_step;}

        xvals[i] = h;
        Y[0][i] = def(h);
        Y[1][i] = util(1, 2, h);
        Y[2][i] = util(1, 3, h);
        Y[3][i] = util(1, 4, h);
        Y[4][i] = util(1, 5, h);
        
        cout << i << endl;
    }

    string names[5] = {
        "Default.temp", 
        "g12.temp",
        "g13.temp",
        "g14.temp",
        "g15.temp"
    };
    
    generate_plot_mpf(5, xvals, Y, names);
}