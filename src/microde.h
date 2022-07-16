#ifndef MICRO_ODE_TOOLKIT_H
#define MICRO_ODE_TOOLKIT_H

/*
 * microde header file by Bram Rodgers.
 * Original Draft Dated: 14, July 2022
 */

/*
 * Header File Body:
 */

/*
 * Macros and Includes go here.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#define mcrd_int long
#define mcrd_flt double
#define mcrd_eps 2.2204e-16
/*
 * Struct Definitions:
 */
typedef struct{
    mcrd_flt* c;//coordinate array for a vector.
    mcrd_int  n;//number of coordinates.
} mcrd_vec;

/*
 * Function Declarations:
 */

/*
 * x = malloc an mcrd vec.
 * x->c = malloc sizeof numel many mcrd_ints
 * x->n = numel
 * return x
 */
mcrd_vec* mcrd_alloc_vec(mcrd_int numel);
void      mcrd_free_vec(mcrd_vec* x);
void mcrd_copy(mcrd_vec* dest, mcrd_vec* src);


mcrd_vec* mcrd_alloc_block(mcrd_int nvecs, mcrd_int numel);

/*
 * Linear combination of 2 vectors with n entries.
 * z->c[i] = a*x->c[i] + b*y->c[i]
 * For each i = 0,1,...,n-1
 */
void mcrd_axpby(mcrd_flt a, mcrd_vec* x, mcrd_flt b, mcrd_vec* y, mcrd_vec* z);

/*
 * Linear combination of 3 vectors with n entries.
 * w->c[i] = a*x->c[i] + b*y->c[i] + c*z->c[i]
 * For each i = 0,1,...,n-1
 */
void mcrd_axpbypcz(mcrd_flt a, mcrd_vec* x,
                   mcrd_flt b, mcrd_vec* y,
                   mcrd_flt c, mcrd_vec* z,
                   mcrd_vec* w);

/*
 * Linear combination with variable argument list.
 * variable argument list expected in alternating order, floats then vecs.
 * Example:
 *      mcrd_lincombo(x, 3, 3.0,z, 4.0,w, -5.0,p);
 * Executes as:
 *          x->c[i] = 3.0*z->c[i] + 4.0*w->c[i] + (-5.0)*p->c[i]
 * 
 * Here, we have numTerms = 3, a sum with 3 terms.
 * Note: In the example above, the sum is started by setting
 *       x->c[i] = 3.0*z->c[i]
 *       then parsing the remaining arguments in a for loop.
 *       If behavior like:
 *       x = x + y
 *       is desired, then  x must be the first term in the variable argument
 *       list.
 */
void mcrd_lincombo(mcrd_vec* x, mcrd_int numTerms, ...);

/*
 * root mean square difference of x and y. (2 norm difference.
 */
mcrd_flt mcrd_mse(mcrd_vec* x, mcrd_vec* y);

//rms norm of input vector.
mcrd_flt mcrd_rms(mcrd_vec* x);

/*
* Single time step of Euler forward implemented using variable argument
* functions. Used for solving an ODE of the form
* 
*  dx
* ----  = vecField(x)
*  dt
* 
* Computes the order 1 estimate
*
* x_new = x_old + dt*vecField(x_old)
*
*
* x_old:      mcrd vector of the current step.
* x_new:      mcrd vector for order 1 estimate of the next step.
* vFld_old:   stores the value of vecField(x_old).
* vecField:   A function pointer where first argument is input and second
*             argument is the return value for an ODE's velocity field.
*             Third variable is the variable argument count, always assumed
*             zero at runtime. It is assumed that static variables
*             are used within the implementation of vecField. These are
*             set using the variable argument count, then not modified
*             during calls by this routine.
* dt:         Time step size for ODE solver.
*/
void mcrd_euler_step(mcrd_vec* x_old,
                     mcrd_vec* x_new,  
                     mcrd_vec* vFld_old,
                     mcrd_flt  dt,
                     void (*vecField)(mcrd_vec*,mcrd_vec*,int,...));

/*
* Single time step of Euler forward implemented using variable argument
* functions. Used for solving an ODE of the form
* 
*  dx
* ----  = vecField(x)
*  dt
* 
* Computes the order 1 estimate
*
* x_stg = x_old + dt*vecField(x_old)
* x_new = x_old + 0.5*dt*(vecField(x_old) + vecField(x_stg))
*
*
* x_old:      mcrd vector of the current step.
* x_new:      mcrd vector for order 2 estimate of the next step.
* vFld_old:   stores the value of vecField(x_old)
* x_stg:      stores the value of x_stg
* vecField:   A function pointer where first argument is input and second
*             argument is the return value for an ODE's velocity field.
*             Third variable is the variable argument count, always assumed
*             zero at runtime. It is assumed that static variables
*             are used within the implementation of vecField. These are
*             set using the variable argument count, then not modified
*             during calls by this routine.
* dt:         Time step size for ODE solver.
*/
void mcrd_heun_step( mcrd_vec* x_old,
                     mcrd_vec* x_new,
                     mcrd_vec* vFld_old,
                     mcrd_vec* x_stg,
                     mcrd_flt  dt,
                     void (*vecField)(mcrd_vec*,mcrd_vec*,int,...));

/*
 * Use a 2 stage embedded Runge-Kutta method to automatically select
 * the a time step satisfying the given absolute tolerance.
 */
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
                      mcrd_flt  abs_tol);

/*
 * x_snap[0] assumed to have t_len many mcrd_vecs allocated to it.
 * to read time snapshot k, access x_snap[0][k].c, x_snap[0][k].n
 * workVec assumed to have 5*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o1(mcrd_vec* x_init,
                       mcrd_vec** x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,   
                       mcrd_flt* workVec);
                       
                       
/*
 * Use a 3 stage embedded Runge-Kutta method to automatically select
 * the a time step satisfying the given absolute tolerance.
 */
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
                      mcrd_flt  abs_tol);
/*
 * x_snap[0] assumed to have t_len many mcrd_vecs allocated to it.
 * to read time snapshot k, access x_snap[0][k].c, x_snap[0][k].n
 * workVec assumed to have 7*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o2(mcrd_vec* x_init,
                       mcrd_vec** x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,
                       mcrd_flt* workVec);

/*
 * Use a 3 stage embedded Runge-Kutta method to automatically select
 * the a time step satisfying the given absolute tolerance.
 */
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
                      mcrd_flt  abs_tol);

/*
 * x_snap[0] assumed to have t_len many mcrd_vecs allocated to it.
 * to read time snapshot k, access x_snap[0][k].c, x_snap[0][k].n
 * workVec assumed to have 9*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o4(mcrd_vec* x_init,
                       mcrd_vec** x_snap,
                       mcrd_flt* t,
                       mcrd_int  t_len,
                       void (*vecField)(mcrd_vec*,mcrd_vec*,int,...),
                       mcrd_flt  abs_tol,
                       mcrd_flt* workVec);

#endif
