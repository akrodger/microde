/*
 * C Source file Template by Bram Rodgers.
 * Original Draft Dated: 25, Feb 2018
 */

/*
 * Macros and Includes go here: (Some common ones included)
 */
#include "microde.h"
#define mcrd_min(a,b) ((a) < (b)) ? (a) : (b)
#define mcrd_max(a,b) ((a) > (b)) ? (a) : (b)
/*
 * Locally used helper functions:
 */

//A few lines taken from ODEs II stiff problems pg 35 by Hairer and Wanner
//for the purposes of selecting the time step of an ode solver.
mcrd_flt mcrd_step_select(mcrd_flt dt,
                          mcrd_flt err_old,
                          mcrd_flt err_new,
                          mcrd_flt abs_tol){
    mcrd_flt step_select;
    step_select = 0.9*pow(abs_tol/err_new,0.35)*pow(err_old/abs_tol,0.2);
    step_select = mcrd_max(0.25,step_select);
    step_select = mcrd_min(1.5,step_select);
    step_select *= dt;
    return step_select;
}
/*
 * Static Local Variables:
 */

/*
 * Function Implementations:
 */

mcrd_vec* mcrd_alloc_vec(mcrd_int numel){
    mcrd_vec* x = (mcrd_vec*) malloc(sizeof(mcrd_vec));
    x->c = (mcrd_flt*) malloc(sizeof(mcrd_flt)*numel);
    x->n = numel;
    return x;
}

mcrd_vec* mcrd_alloc_block(mcrd_int nvecs, mcrd_int numel){
    mcrd_int k = 0;
    mcrd_vec* vec_block;
    vec_block = (mcrd_vec*) malloc(sizeof(mcrd_vec)*nvecs);
    for(k=0;k<nvecs;k++){
        vec_block[k].c = (mcrd_flt*) malloc(sizeof(mcrd_flt)*numel);
        vec_block[k].n = numel;
    }
    return vec_block;
}

void mcrd_free_vec(mcrd_vec* x){
    free(x->c);
    free(x);
}

void mcrd_copy(mcrd_vec* dest, mcrd_vec* src){
    mcrd_int i;
    if(src->n == dest->n){
        for(i=0;i<src->n;i++){
            dest->c[i] = src->c[i];
        }
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_copy::"
                        "dest->n = %ld, src->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) dest->n, (long) src->n);
    }
}


void mcrd_axpby(mcrd_flt a, mcrd_vec* x, mcrd_flt b, mcrd_vec* y, mcrd_vec* z){
    mcrd_int i;
    if(x->n == y->n && y->n == z->n){
        for(i=0;i<x->n;i++){
            z->c[i] = a*x->c[i]+b*y->c[i];
        }
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_axpby::"
                        "x->n = %ld, y->n = %ld, z->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n, (long) z->n);
    }
}

void mcrd_axpbypcz(mcrd_flt a, mcrd_vec* x,
                   mcrd_flt b, mcrd_vec* y,
                   mcrd_flt c, mcrd_vec* z,
                   mcrd_vec* w){
    mcrd_int i;
    if(x->n == y->n && y->n == z->n && z->n == w->n){
        for(i=0;i<x->n;i++){
            w->c[i] = a*x->c[i]+b*y->c[i]+c*z->c[i];
        }
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_axpbypcz::"
                        "x->n = %ld, y->n = %ld, z->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n, (long) z->n);
    }
    
}

mcrd_flt mcrd_mse(mcrd_vec* x, mcrd_vec* y){
    mcrd_int i;
    mcrd_flt mse = 0.0;
    mcrd_flt diff;
    if(x->n != y->n){
        fprintf(stderr, "%s:Line %d %s::mcrd_mse::"
                        "x->n = %ld, y->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n);
    }
    for(i=0;i<x->n;i++){
        diff = (x->c[i] - y->c[i]);
        diff *= diff;
        mse += diff/x->n;
    }
    return sqrt(mse);
}

//rms norm of input vector.
mcrd_flt mcrd_rms(mcrd_vec* x){
    mcrd_int i;
    mcrd_flt rms = 0.0;
    mcrd_flt sqr;
    for(i=0;i<x->n;i++){
        sqr = x->c[i]*x->c[i];
        rms += sqr/x->n;
    }
    return sqrt(rms);
}

void mcrd_euler_step(mcrd_vec* x_old,
                     mcrd_vec* x_new,  
                     mcrd_vec* vFld_old,  
                     mcrd_flt  dt,
                     void (*vecField)(mcrd_vec*,mcrd_vec*,int,...)){
    vecField(x_old,vFld_old,0);
    mcrd_axpby(1.0,x_old,dt,vFld_old,x_new);
}

void mcrd_heun_step( mcrd_vec* x_old,
                     mcrd_vec* x_new,
                     mcrd_vec* vFld_old,
                     mcrd_vec* x_stg,
                     mcrd_flt  dt,
                     void (*vecField)(mcrd_vec*,mcrd_vec*,int,...)){
    mcrd_int i;
    mcrd_euler_step(x_old,x_stg,vFld_old,dt,vecField);
    vecField(x_stg,x_new,0);
    mcrd_axpbypcz(1.0,x_old,0.5*dt,vFld_old,0.5*dt,x_new,x_new);
}


void mcrd_o1_autostep(mcrd_vec* x_old,
                      mcrd_vec* x_new,
                      mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new,
                      mcrd_vec* x_cmp,
                      mcrd_flt* dt_old_ptr,
                      mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                      mcrd_flt* err_old_ptr,
                      mcrd_flt* err_new_ptr,
                      mcrd_flt  abs_tol){
    mcrd_int i = 0;
    mcrd_flt err_old = err_old_ptr[0];
    mcrd_flt err_new = 1.0;
    mcrd_flt dt = *dt_old_ptr;
    while(err_new > abs_tol){
        mcrd_axpby(1.0, x_old, dt, vFld_old, x_new);
        vecField(x_new,vFld_new,0);
        mcrd_axpbypcz(1.0,x_old,0.5*dt,vFld_old,0.5*dt,vFld_new,x_cmp);
        err_old = err_new;
        err_new = mcrd_mse(x_new, x_cmp);
        *dt_old_ptr = dt;
        //compute new dt via handcrafted artisinal PI feedback.
        dt = mcrd_step_select(dt, err_old, err_new, abs_tol);
    }
    err_old_ptr[0] = err_old;
    err_new_ptr[0] = err_new;
    dt_new_ptr[0] = dt;
}

void mcrd_ode_solve_o1(mcrd_vec* x_init,
                       mcrd_vec* x_old,
                       mcrd_vec* x_new,
                       mcrd_vec* vFld_old,
                       mcrd_vec* vFld_new,
                       mcrd_vec* x_cmp,
                       mcrd_vec** x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol){
    static const mcrd_flt eps1 = 1e-5, eps2 = 1e-6, eps3 = 1e-15, eps4=1e-3;
    mcrd_int k = 0, breakflag = 0;
    mcrd_flt n1, n2, n3, dt1, dt2, dt, err_old, err_new, time_now, time_new;
    mcrd_flt theta = 0.0;
    err_old = 1.0;
    err_new = 1.0;
    //initializing vectors.
    mcrd_copy(&(x_snap[0][0]), x_init);
    mcrd_copy(x_old, x_init);
    vecField(x_init,vFld_old,0);
    //begin 1st time step heuristic from Hairer and Wanner, nontstiff ODEs
    n1 = mcrd_rms(&(x_snap[0][0]));
    n2 = mcrd_rms(vFld_old);
    if(n1 < eps1 || n2 < eps1){
        dt1 = eps2;
    }else{
        dt1 = 0.01*n1/n2;
    }
    //do an euler step with this trial time step
    mcrd_axpby(1.0,x_old,dt1,vFld_old,x_new);
    vecField(x_new,vFld_new,0);
    //use that euler step to estimate the ODE's acceleration
    n3 = mcrd_mse(vFld_old,vFld_new);
    if(mcrd_max(n2,n3) > eps3){
        dt2 = sqrt(0.01/mcrd_max(n2,n3));
    }else{
        dt2 = mcrd_max(eps2,dt1*eps4);
    }
    dt = mcrd_min(100*dt1,dt2)*sqrt(abs_tol);
    dt1 = dt;
    //begin main time stepping loop.
    time_now = 0.0;
    k = 1;
    while(breakflag==0){
        mcrd_o1_autostep(x_old,x_new,
                      vFld_old, vFld_new,x_cmp,
                      &dt1,&dt2,vecField,
                      &err_old,&err_new,abs_tol);
        time_new = time_now + dt1;
        if(time_new > t[t_len-1] && k == t_len-1){
            mcrd_axpby(1.0,x_old,
                       t[t_len-1]-time_now,vFld_old,&(x_snap[0][t_len-1]));
            break;
        }
        //use polynomial interponation for dense output.
        while(time_now <= t[k] && t[k] < time_new){
            theta = (t[k] - time_now)/dt1;
            mcrd_axpby(theta,x_old,1.0-theta,x_new,&(x_snap[0][k]));
            k++;
            if(k >= t_len-1){
                breakflag = 1;
                break;
            }
        }
        time_now = time_new;
        mcrd_copy(x_old, x_new);
        mcrd_copy(vFld_old, vFld_new);
        dt1 = dt2;
    }
}
