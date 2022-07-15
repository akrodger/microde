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

mcrd_flt* linspace(mcrd_flt left, mcrd_flt right, mcrd_int n);
/*
 * Template main:
 */
int main(int argc, char** argv){
    char* pEnd = NULL; 
    mcrd_flt pi      = acos(-1.0);
    mcrd_flt t_final = 2*pi;
    mcrd_flt abs_tol = 1e-6;
    mcrd_int Nt;
    if(argc > 1){
        Nt = (mcrd_int) strtol(argv[1], &pEnd, 0);
    }else{
        Nt = 100;
    }
    mcrd_int k       = 0;
    mcrd_flt a       = 1.0;
    mcrd_flt dt      = t_final/Nt;
    mcrd_flt* t       = linspace(0.0, t_final, Nt+1);
    mcrd_vec* x_init   = mcrd_alloc_vec(2);
    mcrd_vec* x_old    = mcrd_alloc_vec(2);
    mcrd_vec* x_new    = mcrd_alloc_vec(2);
    mcrd_vec* vFld_old = mcrd_alloc_vec(2);
    mcrd_vec* vFld_new = mcrd_alloc_vec(2);
    mcrd_vec* x_cmp    = mcrd_alloc_vec(2);
    mcrd_vec* x_snap   = mcrd_alloc_block(Nt+1, 2);
    x_init->c[0] = 1;
    x_init->c[1] = 0;
    mcrd_copy(x_old, x_init);
    vecField(NULL,NULL,1,a);
    mcrd_ode_solve_o1(x_init,x_old,x_new,vFld_old,vFld_new,x_cmp,&x_snap,
                       t, Nt+1,&vecField,  abs_tol);
    for(k=0;k<Nt;k++){
        x_init->c[0] = cos(a*t[k]);
        x_init->c[1] = sin(a*t[k]);
        printf("\nMSE = %le",(double) mcrd_mse(x_init,&(x_snap[k])));
    }
    printf("\n");
    mcrd_free_vec(x_init);
    mcrd_free_vec(x_old);
    mcrd_free_vec(x_new);
    mcrd_free_vec(vFld_old);
    mcrd_free_vec(vFld_new);
    mcrd_free_vec(x_cmp);
    free(t);
    for(k=0;k<=Nt;k++){
        free(x_snap[k].c);
    }
    free(x_snap);
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


mcrd_flt* linspace(mcrd_flt left, mcrd_flt right, mcrd_int n){
    mcrd_int k;
    mcrd_flt* x = (mcrd_flt*) malloc(sizeof(mcrd_flt)*n);
    for(k = 0; k<n; k++){
        x[k] = left + ((k*(right-left))/(n-1));
    }
    return x;
}
