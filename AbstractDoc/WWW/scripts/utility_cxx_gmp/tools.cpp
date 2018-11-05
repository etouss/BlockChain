#include "tools.hpp"

using namespace std;

mpf_class a = 0.9999966993_mpf;
//mpf_class a = 0.9999967030_mpf;
mpf_class b = 0.9999996156_mpf;


// Cat numbers generating function: (1-sqrt(1-4x))/(2x)
mpf_class Cat(mpf_class x){

    mpf_class res;
    res = 1 - 4 * x;

    if(cmp(res,0)<0){
        cout << "Negative sqrt" << endl;
        cout << "Fatal error. Exiting. " << endl;
        exit(0);
    }

    res = 2 / (sqrt(res) + 1);

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


mpz_class fact(int x){
    mpz_class res = 1;
    if (x < 0) {
        cout << "Negative factorial: "<< x << endl;
        exit(0);
        return 0;
    }
    if (x == 0) return 1;
    if (x > 0){
        return x*fact(x-1);
    }
    return -1;
}


// Choose y elements from x
mpz_class choose(int x, int y){
    if(y > x) return 0;
    mpz_class num = fact(x);
    mpz_class den = fact(y) * fact(x-y);

    return num / den;
}


// Domagoj put the name. This is sigma_{a,b}
mpz_class choose_gen(int i, int j){
    if((i < 0) || (j < 0)) return 0;

    mpz_class res = (i + 1) * choose(i + 2 * j, i + j) / (i + j + 1);
    return res;
}

// See better_code.py
mpf_class Cat_g(mpf_class x, int aa, int bb){
    mpf_class sum = 0;
    mpf_class power_of_x = 1;

    for (int i = 0; i < bb + 1; i++){
        sum += power_of_x * (aa + 1) / (aa + i + 1) \
               * choose(aa + 2 * i, aa + i); 
        power_of_x *= x;
    }

    return sum;
}


// See better_code.py. 
mpz_class f(int aa, int r, int bb){
    mpz_class res = 0;
    
    for (int i = 0; i < r - aa + 1; i++)
    {
        res += choose_gen(aa - 1, i) * choose_gen(aa + bb - r, r - aa - i);
    }
    return res;
}

// pentagon generation
mpz_class Pent_aux(int aa, int bb, int r){

    mpz_class res = 0;

    if (r <= aa){
        return choose(bb + r, r);
    }

    if (r <= bb){
        res = choose_gen(bb - r, r);
        for(int i = 1; i < aa + 1; i++){
            res += f(i, r, bb);
        }
        return res;
    }
    res = choose_gen(r - bb, bb);
    for (int i = r - bb + 1; i < aa + 1; i ++){
        res += f(i, r, bb);        
    }

    return res;
}


// real Pent
mpz_class Pent(int aa, int bb, int r){
    return Pent_aux(aa, bb - 1, r);
}


mpf_class a1(int k, int j, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;
    mpf_class res = (a * b * h) / ((1 - a) * (1 - b)) * custom_pow(x, j) * \
        (Cat_g(x, j - 1, l - j) - (custom_pow(a,j + 1) * Cat_g(y, j - 1, l - j)));
    return res;
}


mpf_class sum_a1(int k, int l, mpf_class h){
    mpf_class sum = 0;
    for(int j = 1; j < k + 1; j ++){
        sum += a1(k, j, l, h);
    } 
    return sum;
}


mpf_class a2(int k, int j, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = (a * b * h) * custom_pow(y, j) * Cat_g(y, j - 1, l - j);

    return res;

}


mpf_class sum_a2(int k, int l, mpf_class h){
    mpf_class sum = 0;
    for (int j = 1; j < k + 1; j++){
        sum += a2(k, j, l, h);
    }
    return sum;
}


mpf_class b1(int k, int j, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = custom_pow(b, j + k + 1) * custom_pow(a, j - k + 1) * custom_pow(h, k + 1) * custom_pow(1 - h, j) / ((1 - a) * (1 - b));
    res *= (Cat_g(x, k - 1, l - k) - (custom_pow(a,k + 1) * Cat_g(y, k - 1, l - k)));

    return res;
}


mpf_class sum_b1(int k, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = a * a * custom_pow(x, k + 1) / ((1 - a) * (1 - b) * (1 - a * b * (1 - h))) \
        * (Cat_g(x, k - 1, l - k) - custom_pow(a, k + 1) * Cat_g(y, k - 1, l - k));

    return res;
}
    

mpf_class b2(int k, int j, int l, mpf_class h){
    
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = mpf_class(b, j + k + 1) * custom_pow(a, j + 1) * custom_pow(h, k + 1) * custom_pow(1 - h, j) \
        * Cat_g(y, k - 1, l - k);

    return res;

}


mpf_class sum_b2(int k, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;
    mpf_class res = a * custom_pow(y, k + 1) / (1 - a * b * (1 - h)) * Cat_g(y, k - 1, l - k);

    return res;
}


// E_a,b - Amount of trapezoidal paths
mpf_class E(int aa, int bb){

    mpz_class res = choose(aa + 2 * bb, aa + bb);
    res *= (aa + 1);
    res /= (aa + bb + 1);

    return res;
}
    

mpf_class c(int k, int j, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;
    mpf_class res = 0;

    for(int r = 0; r < l; r++){
        res += custom_pow(b * h, r) * Pent(j - 1, l + 1 - j, r);
    }

    res *= y * custom_pow(a * b * (1 - h), l);

    return res;
}
    

mpf_class sum_c(int k, int l, mpf_class h){
    mpf_class sum = 0;
    for(int j = 1; j < k + 1; j ++){
        sum += c(k, j, l, h);
    } 

    return sum;
}


mpf_class PP(int l, int j, mpf_class x){
    mpf_class res = 0;
    for(int r = 0; r < l; r++){
        res += custom_pow(x, r) * Pent(j - 1, l + 1 - j, r);
    }
    return res;
}



mpf_class d(int k, int j, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = 0;

    for(int r = 0; r < l; r++){
        res += custom_pow(b * h, r) * Pent(k - 1, l + 1 - k, r);
    }
    
    res *= y * custom_pow(a * b * (1 - h), j + l - k);
    return res;
}


mpf_class sum_d(int k, int l, mpf_class h){
    mpf_class x = b * b * h * (1 - h);
    mpf_class y = a * x;

    mpf_class res = y * custom_pow(a * b * (1 - h), l + 1) / (1 - a * b * (1 - h)) \
        * PP(l, k, b * h);

    return res;
}


mpf_class util(int k, int l, mpf_class h){


    mpf_class res = a * b * h / (1 - b) + sum_a1(k, l, h) + sum_b1(k, l, h);

    mpf_class den = 1 - a * b * h - sum_a2(k, l, h) - sum_b2(k, l, h) \
        - sum_c(k, l, h) - sum_d(k, l, h);

    res /= den;

    return (1-b) * res ;
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
    if (j == 0) return a * b * h / (1-b);

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


void generate_plot_mpf(int length, double xvals[NUM_POINTS], mpf_class yvals [][NUM_POINTS], string * names){

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
            fprintf(temp[i], "%lf %lf \n", xvals[j], yvals[i][j].get_d());
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