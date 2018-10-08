#include "btc.hpp"

using namespace std;

string str_a = "0.9999999999999999";
string str_b = "0.99999999999999";

int main(int argc, char** argv){

    mpf_set_default_prec (PREC);
    
    // Bitcoin's alpha:
    mpf_t half;
    mpf_init_set_str(half, "0.5", 10);
    for(int i = 0; i < 17; i++){
        mpf_sqrt(half, half);
    }
    gmp_printf("a = %.*Ff\n", PREC/3, half);



    // Init parameters
    mpf_init(a);
    mpf_init(b);

    // Set parameters from string
    mpf_set_str(a, str_a.c_str(), 10);
    mpf_set_str(b, str_b.c_str(), 10);

    mpf_t h; 
    mpf_init(h);

    mpf_set_str(h, "0.50001", 10);

    mpf_t x;
    mpf_init(x);


    gmp_printf("a = %.*Ff\n", PREC/3, a);
    gmp_printf("b = %.*Ff\n", PREC/3, a);





    cout << "Default :" << endl; 
    def(x, h);

    gmp_printf("%.*Ff\n", PREC/3, x);
    cout << "Always fork :" << endl; 

    always_fork(x, h);
    gmp_printf("%.*Ff\n", PREC/3, x);


    
    //gmp_printf ("fixed point mpf %.*Ff with %d digits\n", n, f, n);

    return 0;
}


void Catalan(mpf_t res, mpf_t x){

    // sqr = sqrt(1-4x)
    mpf_mul_ui(res, x, 4);
    mpf_neg(res, res);
    mpf_add_ui(res, res, 1);
    if(mpf_cmp_ui(res,0)<0){
        cout << "negative sqrt" << endl;
    }
    mpf_sqrt(res, res);

    mpf_neg(res, res);
    mpf_add_ui(res, res, 1);
    mpf_div_ui(res, res, 2);

    mpf_div(res, res, x);

    return;
}

void def(mpf_t res, mpf_t h){

    mpf_t num, den, aux;
    mpf_inits(num, den, aux, NULL);

    mpf_mul(num, a, b);
    mpf_mul(num, num, h);

    mpf_ui_sub(den, 1, b);
    mpf_mul(aux, a, b);
    mpf_neg(aux, aux);
    mpf_add_ui(aux, aux, 1);

    mpf_mul(den, den, aux);

    mpf_div(res, num, den);

    mpf_clears(aux, num, den, NULL);
}


void psi(mpf_t res, mpf_t h){


    mpf_t fact, arg, cat;
    mpf_inits(fact, arg, cat, NULL);

    mpf_mul(fact, a, b);
    mpf_mul(fact, fact, h);

    mpf_ui_sub(arg, 1, h);
    mpf_mul(arg, arg, h);
    mpf_mul(arg, arg, b);
    mpf_mul(arg, arg, b);
    mpf_mul(arg, arg, a);

    Catalan(cat, arg);

    mpf_mul(res, fact, cat);

    mpf_clears(fact, arg, cat, NULL);
}

void phi(mpf_t res, mpf_t h){


    // b2h1mh = b^2 h (1-h)
    mpf_t aux, b2h1mh, fact, num, den, arg, cat1, cat2;
    mpf_inits(aux, b2h1mh, fact, num, den, arg, cat1, cat2, NULL);

    // Auxiliary: 
    mpf_ui_sub(b2h1mh, 1, h);
    mpf_mul(b2h1mh, b2h1mh, h);
    mpf_mul(b2h1mh, b2h1mh, b);
    mpf_mul(b2h1mh, b2h1mh, b);

    // Factor:
    
    mpf_mul(num, a, b);
    mpf_mul(num, num, h);

    mpf_ui_sub(den, 1, a);
    mpf_ui_sub(aux, 1, b);
    mpf_mul(den, den, aux);

    mpf_div(fact, num, den);

    // First Catalan:
    Catalan(cat1, b2h1mh);

    // Second Catalan * a:
    mpf_mul(arg, b2h1mh, a);
    Catalan(cat2, arg);
    mpf_mul(cat2, cat2, a);

    // Difference of Catalans, store in cat1:
    mpf_sub(cat1, cat1, cat2);

    // Result:
    mpf_mul(res, cat1, fact);

    // Free memory
    mpf_clears(aux, b2h1mh, fact, num, den, arg, cat1, cat2, NULL);
}


void always_fork(mpf_t res, mpf_t h){

    mpf_t aux;
    mpf_init(aux);

    psi(aux, h);
    mpf_ui_sub(aux, 1, aux);

    phi(res, h);

    mpf_div(res, res, aux);

    mpf_clear(aux);

}
