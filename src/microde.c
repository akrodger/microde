/*
 * microde source file by Bram Rodgers.
 * Original Draft Dated: 14, July 2022
 */

/*
 * Macros and Includes go here: (Some common ones included)
 */
#include "microde.h"
#if defined(SHMEM_PARA_MICRODE)
#include "omp.h"
#endif
#define mcrd_min(a,b) (((a) < (b)) ? (a) : (b))
#define mcrd_max(a,b) (((a) > (b)) ? (a) : (b))
#define MCRD_FIRST_DT_SETUP \
    {\
    /*begin 1st time step heuristic from Hairer and Wanner, nontstiff ODEs*/\
    n1 = mcrd_rms(&(x_snap[0]));\
    n2 = mcrd_rms(&vFld_old);\
    if(n1 < eps1 || n2 < eps1){\
        dt1 = eps2;\
    }else{\
        dt1 = 0.01*n1/n2;\
    }\
    /*do an euler step with this trial time step*/\
    mcrd_axpby(1.0,&x_old,dt1,&vFld_old,&x_new);\
    vecField(&x_new,&vFld_new,0);\
    /*use that euler step to estimate the ODE's acceleration*/\
    n3 = mcrd_mse(&vFld_old,&vFld_new);\
    if(mcrd_max(n2,n3) > eps3){\
        dt2 = sqrt(0.01/mcrd_max(n2,n3));\
    }else{\
        dt2 = mcrd_max(eps2,dt1*eps4);\
    }\
    dt = mcrd_min(100*dt1,dt2)*sqrt(abs_tol);\
    dt1 = mcrd_max(dt,100*mcrd_eps);\
    }
/*
 * Locally used helper functions:
 */

//helper function for evaluating the bogacki-shampine method.
void mcrd_eval_bgsp(mcrd_vec* x_old,
                    mcrd_vec* x_new,
                    mcrd_vec* x_cmp,
                    mcrd_vec* vFld_old,
                    mcrd_vec* vFld_new,
                    mcrd_vec* vFld_stg1,
                    mcrd_vec* vFld_stg2,
                    mcrd_flt dt,
                    void (*vecField)(mcrd_vec*,mcrd_vec*,int,...)){
    //Butcher tableau constants for Bogacki-Shampine method
    //Butcher matrix:
    static const mcrd_flt bgsp_a_21 = 0.5;
    static const mcrd_flt bgsp_a_32 = 0.75;
    //output weights:
    static const mcrd_flt bgsp_b_1  = 2.0/9.0;
    static const mcrd_flt bgsp_b_2  = 1.0/3.0;
    static const mcrd_flt bgsp_b_3  = 4.0/9.0;
    //output weights for comparison estimate:
    static const mcrd_flt bgsp_e_1  = 7.0/24.0;
    static const mcrd_flt bgsp_e_2  = 0.25;
    static const mcrd_flt bgsp_e_3  = 1.0/3.0;
    static const mcrd_flt bgsp_e_4  = 0.125;
    static mcrd_flt* ptrWrk[5];
    static mcrd_flt scalWrk[5];
    mcrd_vec* x_stg = x_cmp;
    mcrd_axpby(1.0, x_old, dt*bgsp_a_21, vFld_old, x_stg);
    vecField(x_stg,vFld_stg1,0);
    mcrd_axpby(1.0, x_old, dt*bgsp_a_32, vFld_stg1, x_stg);
    vecField(x_stg,vFld_stg2,0);
    mcrd_lincombo(x_new, ptrWrk, scalWrk, 4,      1.0, x_old,
                                            dt*bgsp_b_1, vFld_old,
                                            dt*bgsp_b_2, vFld_stg1,
                                            dt*bgsp_b_3, vFld_stg2);
    vecField(x_new,vFld_new,0);
    mcrd_lincombo(x_cmp, ptrWrk, scalWrk, 5,        1.0, x_old,
                                            dt*bgsp_e_1, vFld_old,
                                            dt*bgsp_e_2, vFld_stg1,
                                            dt*bgsp_e_3, vFld_stg2,
                                            dt*bgsp_e_4, vFld_new);
}

//helper function for evaluating the Dormand * Prince 5(4) method.
void mcrd_eval_dopr(mcrd_vec* x_old,
                    mcrd_vec* x_new,
                    mcrd_vec* x_cmp,
                    mcrd_vec* vFld_old,
                    mcrd_vec* vFld_new,
                    mcrd_vec* vFld_stg1,
                    mcrd_vec* vFld_stg2,
                    mcrd_vec* vFld_stg3,
                    mcrd_vec* vFld_stg4,
                    mcrd_flt dt,
                    void (*vecField)(mcrd_vec*,mcrd_vec*,int,...)){
    //Butcher tableau constants for Bogacki-Shampine method
    //Butcher matrix:
    static const mcrd_flt dopr_a_21 = 1.0/5.0;
    static const mcrd_flt dopr_a_31 = 3.0/40.0;
    static const mcrd_flt dopr_a_32 = 9.0/40.0;
    static const mcrd_flt dopr_a_41 = 44.0/45.0;
    static const mcrd_flt dopr_a_42 = -56.0/15.0;
    static const mcrd_flt dopr_a_43 = 32.0/9.0;
    static const mcrd_flt dopr_a_51 = 19372.0/6561.0;
    static const mcrd_flt dopr_a_52 = -25360.0/2187.0;
    static const mcrd_flt dopr_a_53 = 64448.0/6561.0;
    static const mcrd_flt dopr_a_54 = -212.0/729.0;
    static const mcrd_flt dopr_a_61 = 9017.0/3168.0;
    static const mcrd_flt dopr_a_62 = -355.0/33.0;
    static const mcrd_flt dopr_a_63 = 46732.0/5247.0;
    static const mcrd_flt dopr_a_64 = 49.0/176.0;
    static const mcrd_flt dopr_a_65 = -5103.0/18656.0;
    //output weights:
    static const mcrd_flt dopr_b_1  = 35.0/384.0;
    static const mcrd_flt dopr_b_3  = 500.0/1113.0;
    static const mcrd_flt dopr_b_4  = 125.0/192.0;
    static const mcrd_flt dopr_b_5  = -2187.0/6784.0;
    static const mcrd_flt dopr_b_6  = 11.0/84.0;
    //output weights for comparison estimate:
    static const mcrd_flt dopr_e_1  = 5179.0/57600.0;
    static const mcrd_flt dopr_e_3  = 7571.0/16695.0;
    static const mcrd_flt dopr_e_4  = 393.0/640.0;
    static const mcrd_flt dopr_e_5  = -92097.0/339200.0;
    static const mcrd_flt dopr_e_6  = 187.0/2100.0;
    static const mcrd_flt dopr_e_7  = 1.0/40.0;
    static mcrd_flt* ptrWrk[7];
    static mcrd_flt scalWrk[7];
    mcrd_vec* x_stg = x_cmp;
    mcrd_vec* vFld_stg_extra = x_new;
    mcrd_axpby(1.0, x_old, dt*dopr_a_21, vFld_old, x_stg);
    vecField(x_stg,vFld_stg_extra,0);
    mcrd_lincombo(x_stg, ptrWrk, scalWrk, 3,         1.0, x_old,
                                            dt*dopr_a_31, vFld_old,
                                            dt*dopr_a_32, vFld_stg_extra);
    vecField(x_stg,vFld_stg1,0);
    mcrd_lincombo(x_stg, ptrWrk, scalWrk, 4,         1.0, x_old,
                                            dt*dopr_a_41, vFld_old,
                                            dt*dopr_a_42, vFld_stg_extra,
                                            dt*dopr_a_43, vFld_stg1);
    vecField(x_stg,vFld_stg2,0);
    mcrd_lincombo(x_stg, ptrWrk, scalWrk, 5,         1.0, x_old,
                                            dt*dopr_a_51, vFld_old,
                                            dt*dopr_a_52, vFld_stg_extra,
                                            dt*dopr_a_53, vFld_stg1,
                                            dt*dopr_a_54, vFld_stg2);
    vecField(x_stg,vFld_stg3,0);
    mcrd_lincombo(x_stg, ptrWrk, scalWrk, 6,         1.0, x_old,
                                            dt*dopr_a_61, vFld_old,
                                            dt*dopr_a_62, vFld_stg_extra,
                                            dt*dopr_a_63, vFld_stg1,
                                            dt*dopr_a_64, vFld_stg2,
                                            dt*dopr_a_65, vFld_stg3);
    vecField(x_stg,vFld_stg4,0);
    mcrd_lincombo(x_new, ptrWrk, scalWrk, 6,        1.0, x_old,
                                            dt*dopr_b_1, vFld_old,
                                            dt*dopr_b_3, vFld_stg1,
                                            dt*dopr_b_4, vFld_stg2,
                                            dt*dopr_b_5, vFld_stg3,
                                            dt*dopr_b_6, vFld_stg4);
    vecField(x_new,vFld_new,0);
    mcrd_lincombo(x_cmp, ptrWrk, scalWrk, 7,        1.0, x_old,
                                            dt*dopr_e_7, vFld_new,
                                            dt*dopr_e_1, vFld_old,
                                            dt*dopr_e_3, vFld_stg1,
                                            dt*dopr_e_4, vFld_stg2,
                                            dt*dopr_e_5, vFld_stg3,
                                            dt*dopr_e_6, vFld_stg4);
}

//compute dense outpute polynomial spline for the DOPRI(5,4) method.
void mcrd_dopr_spline(mcrd_vec* spline,
                    mcrd_vec* x_old,
                    mcrd_vec* vFld_old,
                    mcrd_vec* vFld_new,
                    mcrd_vec* vFld_stg1,
                    mcrd_vec* vFld_stg2,
                    mcrd_vec* vFld_stg3,
                    mcrd_vec* vFld_stg4,
                    mcrd_flt dt,
                    mcrd_flt theta){
    static const mcrd_flt dopr_b_1 = 35.0/384.0;
    static const mcrd_flt dopr_b_3 = 500.0/1113.0;
    static const mcrd_flt dopr_b_4 = 125.0/192.0;
    static const mcrd_flt dopr_b_5 = -2187.0/6784.0;
    static const mcrd_flt dopr_b_6 = 11.0/84.0;
    static mcrd_flt* ptrWrk[7];
    static mcrd_flt scalWrk[7];
    mcrd_flt b1 = (theta*theta)*(3.0-(2.0*theta))*dopr_b_1 +
                  theta*((theta-1)*(theta-1)) -
                  (theta*theta)*((theta-1)*(theta-1))*5.0*(
                  2558722523-(31403016*theta)
                  )/11282082423;
    mcrd_flt b3 = (theta*theta)*(3.0-(2.0*theta))*dopr_b_3 +
                  (theta*theta)*((theta-1)*(theta-1))*100.0*(
                  882725551-(15701508*theta)
                  )/32700410799;
    mcrd_flt b4 = (theta*theta)*(3.0-(2.0*theta))*dopr_b_4 -
                  (theta*theta)*((theta-1)*(theta-1))*25.0*(
                  443332067-(31403016*theta)
                  )/1880347072;
    mcrd_flt b5 = (theta*theta)*(3.0-(2.0*theta))*dopr_b_5 +
                  (theta*theta)*((theta-1)*(theta-1))*32805.0*(
                  23143187-(3489224*theta)
                  )/199316789632;
    mcrd_flt b6 = (theta*theta)*(3.0-(2.0*theta))*dopr_b_6 -
                  (theta*theta)*((theta-1)*(theta-1))*55.0*(
                  29972135-(7076736*theta)
                  )/822651844;
    mcrd_flt b7 = (theta*theta)*(theta-1) +
                  (theta*theta)*((theta-1)*(theta-1))*10.0*(
                  7414447.0-(829305*theta)
                  )/29380423.0;
    mcrd_lincombo(spline, ptrWrk, scalWrk,7,  1.0, x_old,
                                            dt*b1, vFld_old,
                                            dt*b7, vFld_new,
                                            dt*b3, vFld_stg1,
                                            dt*b4, vFld_stg2,
                                            dt*b5, vFld_stg3,
                                            dt*b6, vFld_stg4);
}
//A few lines taken from ODEs II stiff problems pg 28 by Hairer and Wanner
//for the purposes of selecting the time step of an ode solver.
mcrd_flt mcrd_step_select(mcrd_flt dt,
                          mcrd_flt err_old,
                          mcrd_flt err_new,
                          mcrd_flt tol,
                          mcrd_int p){
    mcrd_flt step_select;
    step_select = pow(tol/err_new,1.0/(p+1.0))*pow(err_old/tol,0.08);
    step_select = mcrd_max(0.5,step_select);
    step_select = mcrd_min(2.0,step_select);
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
    vec_block[0].c = (mcrd_flt*) malloc(sizeof(mcrd_flt)*numel*nvecs);
    vec_block[0].n = numel;
    for(k=1;k<nvecs;k++){
        vec_block[k].c = &(vec_block[k-1].c[numel]);
        vec_block[k].n = numel;
    }
    return vec_block;
}

void mcrd_free_vec(mcrd_vec* x){
    free(x->c);
    free(x);
}

void mcrd_copy(mcrd_vec* dest, mcrd_vec* src){
    mcrd_int i, N, loopStart, loopEnd, thNum, thTot, chunkSize, chunkRem;
    thNum = 0;
    thTot = 1;
    if(src->n == dest->n){
        #if defined(SHMEM_PARA_MICRODE)
        #pragma omp parallel default(shared)\
                private(i,chunkSize,N,\
                       loopStart,loopEnd,chunkRem,thNum,thTot)
        {
        thTot = omp_get_num_threads();
        thNum = omp_get_thread_num();
        #endif 
        N = dest->n;
        chunkSize = N / thTot;
        chunkRem  = N % thTot;
        if(thNum < chunkRem){
            loopStart = (chunkSize+1)*thNum;
            loopEnd   = loopStart + chunkSize + 1;
        }else{
            loopStart = (chunkSize+1)*chunkRem + chunkSize*(thNum-chunkRem);
            loopEnd   = loopStart + chunkSize;
        }
        for(i=loopStart;i<loopEnd;i++){
            dest->c[i] = src->c[i];
        }
        #if defined(SHMEM_PARA_MICRODE)
        #pragma omp barrier
        }
        #endif
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_copy::"
                        "dest->n = %ld, src->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) dest->n, (long) src->n);
        exit(1);
    }
}


void mcrd_axpby(mcrd_flt a, mcrd_vec* x, mcrd_flt b, mcrd_vec* y, mcrd_vec* z){
    static mcrd_flt* ptrWrk[2];
    static mcrd_flt  scalWrk[2];
    if(x->n == y->n && y->n == z->n){
        if(z == y){
            mcrd_lincombo(z, ptrWrk, scalWrk, 2, b,y,a,x);
        }else{
            mcrd_lincombo(z, ptrWrk, scalWrk, 2, a,x,b,y);
        }
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_axpby::"
                        "x->n = %ld, y->n = %ld, z->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n, (long) z->n);
        exit(1);
    }
}

void mcrd_axpbypcz(mcrd_flt a, mcrd_vec* x,
                   mcrd_flt b, mcrd_vec* y,
                   mcrd_flt c, mcrd_vec* z,
                   mcrd_vec* w){
    static mcrd_flt* ptrWrk[3];
    static mcrd_flt  scalWrk[3];
    if(x->n == y->n && y->n == z->n && z->n == w->n){
        if(w == z){
            mcrd_lincombo(w, ptrWrk, scalWrk, 3, c,z,a,x,a,y);
        }else if(w == y){
            mcrd_lincombo(w, ptrWrk, scalWrk, 3, b,y,a,x,c,z);
        }else{
            mcrd_lincombo(w, ptrWrk, scalWrk, 3, a,x,b,y,c,z);
        }
    }else{
        fprintf(stderr, "%s:Line %d %s::mcrd_axpbypcz::"
                        "x->n = %ld, y->n = %ld, z->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n, (long) z->n);
        exit(1);
    }
    
}

void mcrd_lincombo(mcrd_vec* x, mcrd_flt** ptrWrk, mcrd_flt* scalWrk,
                   mcrd_int numTerms, ...){
    va_list arglist;
    mcrd_int i,j, N, loopStart, loopEnd, thNum, thTot, chunkSize, chunkRem;
    mcrd_vec* z = NULL;
    mcrd_flt  a = 0.0;
    mcrd_flt  e = 0.0;
    mcrd_flt  y;
    mcrd_flt  t;
    thNum = 0;
    thTot = 1;
    va_start(arglist, numTerms);;
    for(j=0;j<numTerms;j++){
        a = va_arg(arglist, mcrd_flt);
        z = va_arg(arglist, mcrd_vec*);
        scalWrk[j] = a;
        ptrWrk[j] = z->c;
        if(x->n != z->n){
            fprintf(stderr, "%s:Line %d %s::mcrd_lincombo::"
                        "x->n = %ld, term = %ld, z->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) j, (long) z->n);
        exit(1);
        }
    }
    #if defined(SHMEM_PARA_MICRODE)
    #pragma omp parallel default(shared) private(i,j,a,e,y,t,N,\
                       loopStart,loopEnd,chunkRem,thNum,thTot)
    {
    thTot = omp_get_num_threads();
    thNum = omp_get_thread_num();
    #endif 
    N = x->n;
    chunkSize = N / thTot;
    chunkRem  = N % thTot;
    if(thNum < chunkRem){
        loopStart = (chunkSize+1)*thNum;
        loopEnd   = loopStart + chunkSize + 1;
    }else{
        loopStart = (chunkSize+1)*chunkRem + chunkSize*(thNum-chunkRem);
        loopEnd   = loopStart + chunkSize;
    }
    for(i=loopStart;i<loopEnd;i++){
        x->c[i] = 0.0;
        e       = 0.0;
        for(j=0;j<numTerms;j++){
            a       = scalWrk[j];
            t       = x->c[i];
            y       = (a*ptrWrk[j][i]) + e;
            x->c[i] = t + y;
            e       = (t - x->c[i]) + y;
        }
    }
    #if defined(SHMEM_PARA_MICRODE)
    }
    #endif
}

mcrd_flt mcrd_mse(mcrd_vec* x, mcrd_vec* y){
    mcrd_int i, N;
    #if defined(SHMEM_PARA_MICRODE)
	mcrd_int chunkSize;
	#endif
    mcrd_flt mse;
    mcrd_flt diff;
    mse = 0.0;
    if(x->n != y->n){
        fprintf(stderr, "%s:Line %d %s::mcrd_mse::"
                        "x->n = %ld, y->n = %ld",
                 __FILE__, __LINE__, __func__,
                 (long) x->n, (long) y->n);
        exit(1);
    }
    #if defined(SHMEM_PARA_MICRODE)
    #pragma omp parallel default(shared) private(i,diff,chunkSize,N)
    {
    chunkSize = ((x->n) / ((mcrd_int) omp_get_num_threads()))+1;
    #endif        
    N = x->n;
    #if defined(SHMEM_PARA_MICRODE)
    #pragma omp for schedule(dynamic,chunkSize) reduction(+:mse)
    #endif
    for(i=0;i<x->n;i++){
        diff = (x->c[i] - y->c[i]);
        diff *= diff/N;
        mse = mse + diff;
    }
    #if defined(SHMEM_PARA_MICRODE)
    }
    #endif
    return sqrt(mse);
}

//rms norm of input vector.
mcrd_flt mcrd_rms(mcrd_vec* x){
    mcrd_int i, N;
    #if defined(SHMEM_PARA_MICRODE)
	mcrd_int chunkSize;
	#endif
    mcrd_flt rms;
    mcrd_flt sqr;
    rms = 0.0;
    #if defined(SHMEM_PARA_MICRODE)
    #pragma omp parallel default(shared) private(i,sqr,chunkSize,N)
    {
    chunkSize = ((x->n) / ((mcrd_int) omp_get_num_threads()))+1;
    #endif
    N = x->n;
    #if defined(SHMEM_PARA_MICRODE)
    #pragma omp for schedule(dynamic,chunkSize) reduction(+:rms)
    #endif
    for(i=0;i<x->n;i++){
        sqr = x->c[i]*x->c[i]/N;
        rms = rms + sqr;
    }
    #if defined(SHMEM_PARA_MICRODE)
    }
    #endif
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
                      mcrd_flt* nrm_old_ptr,
                      mcrd_flt* nrm_new_ptr,
                      mcrd_flt  abs_tol,
                      mcrd_flt  rel_tol){
    mcrd_flt err_old = err_old_ptr[0];
    mcrd_flt err_new = err_old;
    mcrd_flt nrm_old = nrm_old_ptr[0];
    mcrd_flt nrm_new = 1.0;
    mcrd_flt tol_bound;
    mcrd_flt dt = *dt_old_ptr;
    do{
        mcrd_axpby(1.0, x_old, dt, vFld_old, x_new);
        vecField(x_new,vFld_new,0);
        mcrd_axpbypcz(1.0,x_old,0.5*dt,vFld_old,0.5*dt,vFld_new,x_cmp);
        err_old = err_new;
        err_new = mcrd_mse(x_new, x_cmp);
        nrm_new = mcrd_rms(x_new);
        *dt_old_ptr = dt;
        //compute new dt via handcrafted artisinal PI feedback.
        dt = mcrd_step_select(dt, err_old, err_new, abs_tol,1);
        tol_bound = abs_tol + mcrd_max(nrm_old,nrm_new)*(rel_tol);
    }while(err_new > tol_bound);
    nrm_new_ptr[0] = nrm_new;
    err_old_ptr[0] = err_old;
    err_new_ptr[0] = err_new;
    dt_new_ptr[0] = dt;
}

void mcrd_ode_solve_o1(mcrd_vec* x_init,
                       mcrd_vec* x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,
                       mcrd_flt  rel_tol,
                       mcrd_flt* workVec){
    static const mcrd_flt eps1 = 1e-5, eps2 = 1e-6, eps3 = 1e-15, eps4=1e-3;
    static mcrd_vec x_old, x_new, vFld_old, vFld_new, x_cmp;
    mcrd_int k = 0, breakflag = 0;
    mcrd_flt n1, n2, n3, dt1, dt2, dt, err_old, err_new, time_now, time_new;
    mcrd_flt theta = 0.0;
    x_old.c    = workVec;
    x_old.n    = x_init->n;
    x_new.c    = &(workVec[x_init->n]);
    x_new.n    = x_init->n;
    vFld_old.c = &(workVec[2*x_init->n]);
    vFld_old.n = x_init->n;
    vFld_new.c = &(workVec[3*x_init->n]);
    vFld_new.n = x_init->n;
    x_cmp.c    = &(workVec[4*x_init->n]);
    x_cmp.n    = x_init->n;
    err_old = 1.0;
    err_new = 1.0;
    /*initializing vectors.*/
    mcrd_copy(&(x_snap[0]), x_init);
    mcrd_copy(&x_old, x_init);
    vecField(x_init,&vFld_old,0);
    MCRD_FIRST_DT_SETUP
    //begin main time stepping loop.
    time_now = 0.0;
    k = 1;
    while(breakflag==0){
        mcrd_o1_autostep(&x_old,&x_new,
                      &vFld_old, &vFld_new,&x_cmp,
                      &dt1,&dt2,vecField,
                      &err_old,&err_new,&n1,&n2,abs_tol,rel_tol);
        n1 = n2;
        time_new = time_now + dt1;
        if(time_new > t[t_len-1] && k == t_len-1){
            mcrd_axpby(1.0,&x_old,
                       t[t_len-1]-time_now,&vFld_old,&(x_snap[t_len-1]));
            break;
        }
        //use polynomial interponation for dense output.
        while(time_now <= t[k] && t[k] < time_new){
            theta = (t[k] - time_now)/dt1;
            mcrd_axpby(theta,&x_old,1.0-theta,&x_new,&(x_snap[k]));
            k++;
            if(k >= t_len){
                breakflag = 1;
                break;
            }
        }
        time_now = time_new;
        mcrd_copy(&x_old, &x_new);
        mcrd_copy(&vFld_old, &vFld_new);
        dt1 = dt2;
    }
}

void mcrd_o2_autostep(mcrd_vec* x_old,
                      mcrd_vec* x_new,
                      mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new,
                      mcrd_vec* vFld_stg1,
                      mcrd_vec* vFld_stg2,
                      mcrd_vec* x_cmp,
                      mcrd_flt* dt_old_ptr,
                      mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                      mcrd_flt* err_old_ptr,
                      mcrd_flt* err_new_ptr,
                      mcrd_flt* nrm_old_ptr,
                      mcrd_flt* nrm_new_ptr,
                      mcrd_flt  abs_tol,
                      mcrd_flt  rel_tol){
    mcrd_flt err_old = err_old_ptr[0];
    mcrd_flt err_new = err_old;
    mcrd_flt nrm_old = nrm_old_ptr[0];
    mcrd_flt nrm_new = 1.0;
    mcrd_flt tol_bound;
    mcrd_flt dt = *dt_old_ptr;
    do{
       mcrd_eval_bgsp(   x_old,   x_new,    x_cmp,
                      vFld_old,vFld_new,vFld_stg1,vFld_stg2,
                      dt,vecField);
        err_old = err_new;
        err_new = mcrd_mse(x_new, x_cmp);
        nrm_new = mcrd_rms(x_new);
        *dt_old_ptr = dt;
        //compute new dt via handcrafted artisinal PI feedback.
        dt = mcrd_step_select(dt, err_old, err_new, abs_tol,2);
        tol_bound = abs_tol + mcrd_max(nrm_old,nrm_new)*(rel_tol);
    }while(err_new > tol_bound);
    nrm_new_ptr[0] = nrm_new;
    err_old_ptr[0] = err_old;
    err_new_ptr[0] = err_new;
    dt_new_ptr[0] = dt;
}


void mcrd_ode_solve_o2(mcrd_vec* x_init,
                       mcrd_vec* x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,
                       mcrd_flt  rel_tol,
                       mcrd_flt* workVec){
    static const mcrd_flt eps1 = 1e-5, eps2 = 1e-6, eps3 = 1e-15, eps4=1e-3;
    static mcrd_vec x_old,     x_new, vFld_old, vFld_new,
                    x_cmp, vFld_stg1,vFld_stg2;
    static mcrd_flt* ptrWrk[4];
    static mcrd_flt  scalWrk[4];
    mcrd_int k = 0, breakflag = 0;
    mcrd_flt n1, n2, n3, m1, m2, m3, m4;
    mcrd_flt dt1, dt2, dt, err_old, err_new, time_now, time_new;
    mcrd_flt theta = 0.0;
    x_old.c     = workVec;
    x_old.n     = x_init->n;
    x_new.c     = &(workVec[x_init->n]);
    x_new.n     = x_init->n;
    vFld_old.c  = &(workVec[2*x_init->n]);
    vFld_old.n  = x_init->n;
    vFld_new.c  = &(workVec[3*x_init->n]);
    vFld_new.n  = x_init->n;
    x_cmp.c     = &(workVec[4*x_init->n]);
    x_cmp.n     = x_init->n;
    vFld_stg1.c = &(workVec[5*x_init->n]);
    vFld_stg1.n = x_init->n;
    vFld_stg2.c = &(workVec[6*x_init->n]);
    vFld_stg2.n = x_init->n;
    err_old = 1.0;
    err_new = 1.0;
    /*initializing vectors.*/
    mcrd_copy(&(x_snap[0]), x_init);
    mcrd_copy(&x_old, x_init);
    vecField(x_init,&vFld_old,0);
    MCRD_FIRST_DT_SETUP
    //begin main time stepping loop.
    time_now = 0.0;
    k = 1;
    while(breakflag==0){
        mcrd_o2_autostep(&x_old,&x_new,
                      &vFld_old, &vFld_new,
                      &vFld_stg1, &vFld_stg2,
                      &x_cmp,
                      &dt1,&dt2,vecField,
                      &err_old,&err_new,&n1,&n2,abs_tol,rel_tol);
        n1 = n2;
        time_new = time_now + dt1;
        //use polynomial interponation for dense output.
        while(time_now <= t[k] && t[k] < time_new){
            theta = (t[k] - time_now)/dt1;
            m1 = 1.0-theta-(theta*(theta-1.0)*(1.0-(2.0*theta)));
            m2 =     theta+(theta*(theta-1.0)*(1.0-(2.0*theta)));
            m3 = dt1*theta*(theta-1.0)*(theta-1.0);
            m4 = dt1*theta*(theta-1.0)*theta;
            mcrd_lincombo(&(x_snap[k]),ptrWrk,scalWrk,4,
                           m1,&x_old,m2,&x_new,m3,&vFld_old,m4,&vFld_new);
            k++;
            if(k >= t_len){
                breakflag = 1;
                break;
            }
        }
        time_now = time_new;
        mcrd_copy(&x_old, &x_new);
        mcrd_copy(&vFld_old, &vFld_new);
        dt1 = dt2;
    }
}


void mcrd_o4_autostep(mcrd_vec* x_old,
                      mcrd_vec* x_new,
                      mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new,
                      mcrd_vec* vFld_stg1,
                      mcrd_vec* vFld_stg2,
                      mcrd_vec* vFld_stg3,
                      mcrd_vec* vFld_stg4,
                      mcrd_vec* x_cmp,
                      mcrd_flt* dt_old_ptr,
                      mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                      mcrd_flt* err_old_ptr,
                      mcrd_flt* err_new_ptr,
                      mcrd_flt* nrm_old_ptr,
                      mcrd_flt* nrm_new_ptr,
                      mcrd_flt  abs_tol,
                      mcrd_flt  rel_tol){
    mcrd_flt err_old = err_old_ptr[0];
    mcrd_flt err_new = err_old;
    mcrd_flt nrm_old = nrm_old_ptr[0];
    mcrd_flt nrm_new = 1.0;
    mcrd_flt tol_bound;
    mcrd_flt dt = *dt_old_ptr;
    do{
        mcrd_eval_dopr(   x_old,   x_new,    x_cmp,
                      vFld_old,vFld_new,vFld_stg1,vFld_stg2,
                      vFld_stg3,vFld_stg4, dt,vecField);
        err_old = err_new;
        err_new = mcrd_mse(x_new, x_cmp);
        nrm_new = mcrd_rms(x_new);
        *dt_old_ptr = dt;
        //compute new dt via handcrafted artisinal PI feedback.
        dt = mcrd_step_select(dt, err_old, err_new, abs_tol,4);
        tol_bound = abs_tol + mcrd_max(nrm_old,nrm_new)*(rel_tol);
    }while(err_new > tol_bound);
    nrm_new_ptr[0] = nrm_new;
    err_old_ptr[0] = err_old;
    err_new_ptr[0] = err_new;
    dt_new_ptr[0] = dt;
}


void mcrd_ode_solve_o4(mcrd_vec* x_init,
                       mcrd_vec* x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,
                       mcrd_flt  rel_tol,
                       mcrd_flt* workVec){
    static const mcrd_flt eps1 = 1e-5, eps2 = 1e-6, eps3 = 1e-15, eps4=1e-3;
    static mcrd_vec x_old,     x_new, vFld_old, vFld_new,
                    x_cmp, vFld_stg1, vFld_stg2, vFld_stg3, vFld_stg4;
    mcrd_int k = 0, breakflag = 0;
    mcrd_flt n1, n2, n3;
    mcrd_flt dt1, dt2, dt, err_old, err_new, time_now, time_new;
    mcrd_flt theta = 0.0;
    x_old.c     = workVec;
    x_old.n     = x_init->n;
    x_new.c     = &(workVec[x_init->n]);
    x_new.n     = x_init->n;
    vFld_old.c  = &(workVec[2*x_init->n]);
    vFld_old.n  = x_init->n;
    vFld_new.c  = &(workVec[3*x_init->n]);
    vFld_new.n  = x_init->n;
    x_cmp.c     = &(workVec[4*x_init->n]);
    x_cmp.n     = x_init->n;
    vFld_stg1.c = &(workVec[5*x_init->n]);
    vFld_stg1.n = x_init->n;
    vFld_stg2.c = &(workVec[6*x_init->n]);
    vFld_stg2.n = x_init->n;
    vFld_stg3.c = &(workVec[7*x_init->n]);
    vFld_stg3.n = x_init->n;
    vFld_stg4.c = &(workVec[8*x_init->n]);
    vFld_stg4.n = x_init->n;
    err_old = 1.0;
    err_new = 1.0;
    /*initializing vectors.*/
    mcrd_copy(&(x_snap[0]), x_init);
    mcrd_copy(&x_old, x_init);
    vecField(x_init,&vFld_old,0);
    MCRD_FIRST_DT_SETUP
    //begin main time stepping loop.
    time_now = 0.0;
    k = 1;
    while(breakflag==0){
        mcrd_o4_autostep(&x_old,&x_new,
                      &vFld_old, &vFld_new,
                      &vFld_stg1, &vFld_stg2,
                      &vFld_stg3, &vFld_stg4,
                      &x_cmp,
                      &dt1,&dt2,vecField,
                      &err_old,&err_new,&n1,&n2,abs_tol,rel_tol);
        n1 = n2;
        time_new = time_now + dt1;
        //use polynomial interponation for dense output.
        while(time_now <= t[k] && t[k] < time_new){
            theta = (t[k] - time_now)/dt1;
            mcrd_dopr_spline(&(x_snap[k]),
                             &x_old,
                             &vFld_old,
                             &vFld_new,
                             &vFld_stg1,
                             &vFld_stg2,
                             &vFld_stg3,
                             &vFld_stg4,
                             dt1,
                             theta);
            k++;
            if(k >= t_len){
                breakflag = 1;
                break;
            }
        }
        time_now = time_new;
        mcrd_copy(&x_old, &x_new);
        mcrd_copy(&vFld_old, &vFld_new);
        dt1 = dt2;
    }
}
