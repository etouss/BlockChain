#include "tools.hpp"

using namespace std;


mpf_class a = 0.9997_mpf;
mpf_class b = 0.999997_mpf;


// Cat numbers generating function: (1-sqrt(1-4x))/(2x)
mpf_class Cat(mpf_class x){

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
    cat = Cat(arg);
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

    // First Cat:
    b2h1mh = b * b * h *(1 - h);
    cat1 = Cat(b2h1mh);

    // Second Cat * a:
    arg = a * b2h1mh;
    cat2 = Cat(arg);

    // Result:
    res = (cat1 - a * cat2) * fact ;

    return res;
}


// Auxiliary function for always_fork_k()
mpf_class phi_j(int k, int j, mpf_class h){

    //phi_0 = abh
    if (j == 0) return a * b * h;

    mpf_class x, y, fact, fact2, fact3;
    // x = b^2h(1-h)
    // y = ab^2h(1-h)
    x = b * b * h * (1 - h);
    y = a * x;

    // Fact = abh / (1-a)(1-b)
    fact = a * b * h / ((1 - a) * (1 - b));

    if(j <= k){
        // fact2 = x^j
        fact2 = custom_pow(x,j);
        
        // fact3 = cat(x)^j - a^{j+1} cat(y)^j
        fact3 = custom_pow(Cat(x), j);
        fact3 -= custom_pow(a, j+1) * custom_pow(Cat(y),j);

        return fact * fact2 * fact3;
    }

    // If j > k:

    // fact2 = a^{j-k}b^{j+k}h^k(1-h)^j
    fact2 = 
        custom_pow(a, j-k) *
        custom_pow(b, j+k) *
        custom_pow(h, k)   *
        custom_pow(1-h, j);

    // fact3 = cat(x)^k - a^{k+1} cat(y)^k
    fact3 = custom_pow(Cat(x), k);
    fact3 -= custom_pow(a, k+1) * custom_pow(Cat(y),k);

    return fact * fact2 * fact3;
}


// Auxiliary function for always_fork_k()
mpf_class psi_j(int k, int j, mpf_class h){

    //phi_0 = abh
    if (j == 0) return a * b * h;

    mpf_class x, y, fact, fact2, fact3;
    // x = b^2h(1-h)
    // y = ab^2h(1-h)
    x = b * b * h * (1 - h);
    y = a * x;

    // Fact = abh
    fact = a * b * h;


    if (j <= k){
        // fact2 = y^j
        fact2 = custom_pow(y, j);
        // fact3 = Cat_{j-1}(y)
        fact3 = custom_pow(Cat(y), j);
        return fact * fact2 * fact3;
    }

    // If j > k
    // fact2 = a^{j}b^{j+k}h^k(1-h)^j
    fact2 = 
        custom_pow(a, j)     *
        custom_pow(b, j + k) *
        custom_pow(h, k)     *
        custom_pow(1-h , j);
    // fact3 = Cat_{k-1}(y)
    fact3 = custom_pow(Cat(y), k);
    return fact * fact2 * fact3;

}


mpf_class PHI_k(int k, mpf_class h){
    mpf_class res = 0_mpf;
    for(int j = 0; j < SUM_LIMIT; j++){
        res += phi_j(k, j ,h);
    }
    return res;
}


mpf_class PSI_k(int k, mpf_class h){
    mpf_class res = 0_mpf;
    for(int j = 0; j < SUM_LIMIT; j++){
        res += psi_j(k, j ,h);
    }
    return res;
}


mpf_class PHI_k_closed(int k, mpf_class h){
    
    mpf_class x, y, c_x, c_y, res;
    
    // x = b^2h(1-h)
    // y = ab^2h(1-h)
    x = b * b * h * (1 - h);
    y = a * x;

    c_x = x * Cat(x);
    c_y = y * Cat(y);

    res = 
        (c_x - custom_pow(c_x,k+1)) / (1 - c_x) 
        - a * (c_y - custom_pow(c_y, k+1)) / (1 - c_y)
        + custom_pow(x, k) *a*b*(1-h)* (custom_pow(Cat(x),k)- custom_pow(Cat(y),k)*a)  / (1-a*b*(1-h));
    
    res *= a * b * h / ((1 - a)*(1 - b));

    res += a*b*h/(1-b);

    return res;
}


mpf_class PSI_k_closed(int k, mpf_class h){
    mpf_class x, y, c_y, res;
    
    // x = b^2h(1-h)
    // y = ab^2h(1-h)
    x = b * b * h * (1 - h);
    y = a * x;
    c_y = y * Cat(y);


    res = 
          (c_y - custom_pow(c_y,k + 1))/(1 - c_y) 
        + (a * b * (1-h)) * (custom_pow(c_y,k)) / (1 - a*b*(1-h));

    res *= a * b * h;    
    res += a * b * h;    

    return res;

}


mpf_class window_fork(int k, mpf_class h){
    mpf_class num = PHI_k_closed(k, h);
    mpf_class den = 1 - PSI_k_closed(k, h);
    return (1-b)*num / den;
}




mpf_class custom_pow(mpf_class base, int exp){
    if(exp < 0){
        cout << "Fatal error in " << __func__ << endl;
    }
    mpf_class res;
    mpf_pow_ui(res.get_mpf_t(),base.get_mpf_t(), exp);
    return res;
}


void generate_plot(int length, double xvals[NUM_POINTS], double yvals [][NUM_POINTS], string * names){

    string command = "plot";
    for(int i = 0; i < length; i++){
        command += " '" + names[i] + "' with lines,"; 
    }

    cout << command << endl;

    string commandsForGnuplot[] = {
        "set terminal wxt size 1500,800",
        "unset title", 
        "set key bottom ", 
        "set grid xtics lt 0 lw 1 ",
        "set grid ytics lt 0 lw 1 ",
        command.c_str()
    };
    FILE * temp[length];
    for(int i = 0; i < length; i++){
        temp[i] = fopen((const char *) names[i].c_str(), "w");
        for (int j = 0; j < NUM_POINTS; j++)
        {
            fprintf(temp[i], "%lf %lf \n", xvals[j], yvals[i][j]);
        }
    }
    
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

    for (int i = 0; i < NUM_COMMANDS; i++)
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