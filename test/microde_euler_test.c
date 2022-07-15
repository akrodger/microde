/*
 * C Entry Point file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Macros and Includes go here: (Some common ones listed)
 */
#include <stdio.h>
#include <stdlib.h>
#include "../src/microde.h"
/*
 * Function declarations go here:
 */

void vecField(mcrd_vec* x,mcrd_vec* dxdt,int argc,...);

/*
 * Template main:
 */
int main(int argc, char** argv){
    char* pEnd = NULL; 
    mcrd_flt pi      = acos(-1.0);
    mcrd_flt t_final = pi/2;
    mcrd_int Nt;
    if(argc > 1){
        Nt = (mcrd_int) strtol(argv[1], &pEnd, 0);
    }else{
        Nt = 1000;
    }
    mcrd_int k       = 0;
    mcrd_flt a       = 1.0;
    mcrd_flt dt      = t_final/Nt;
    mcrd_vec* x_init  = mcrd_alloc_vec(2);
    mcrd_vec* x_old   = mcrd_alloc_vec(2);
    mcrd_vec* x_new   = mcrd_alloc_vec(2);
    mcrd_vec* wk      = mcrd_alloc_vec(2);
    x_init->c[0] = 1;
    x_init->c[1] = 0;
    mcrd_copy(x_old, x_init);
    vecField(NULL,NULL,1,a);
    for(k=0;k<Nt;k++){
        mcrd_euler_step(x_old,x_new,wk,dt,&vecField);
        mcrd_copy(x_old, x_new);
    }
    x_init->c[0] = cos(a*t_final);
    x_init->c[1] = sin(a*t_final);
    printf("\nMSE = %le\n",(double) mcrd_mse(x_init,x_new));
    mcrd_free_vec(x_init);
    mcrd_free_vec(x_old);
    mcrd_free_vec(x_new);
    mcrd_free_vec(wk);
    return 0;
}

void vecField(mcrd_vec* x,mcrd_vec* dxdt,int argc,...){
    static mcrd_flt a;
    va_list arglist;
    if(argc == 0){
        dxdt->c[0] =  -a*x->c[1];
        dxdt->c[1] =  a*x->c[0];
    }else if(argc == 1){
        va_start(arglist, argc);
        a = va_arg(arglist, mcrd_flt);
    }
}
