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
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Macro for changing the integer size in use. Default is long.
#define mcrd_int long
// Macro for changing the floating point precision. Default is double.
#define mcrd_flt double
// Change this to be the floating point machine epsilon based on mcrd_float's
// def
#define mcrd_eps 2.2204e-16
// The following macro is used to turn shared memory parallelism using
// OpenMP on and off. An #ifdef check is done on the rms, mse, various
// linear combinaion functions, and copy function.
//#define SHMEM_PARA_MICRODE
/*
 * Struct Definitions:
 */
// A struct for storing the solution to an ODE at a single time step.
typedef struct {
  mcrd_flt* c;  // coordinate array for a vector.
  mcrd_int n;   // number of coordinates.
} mcrd_vec;

/*
 * Function Declarations:
 */

/*
 * Memory allocation for the mcrd_vec struct.
 * Blueprint:
 *  x = malloc an mcrd vec.
 *  x->c = malloc sizeof numel many mcrd_ints
 *  x->n = numel
 *  return x
 *
 * Args:
 *  numel:    Number of elemnts for vector.
 *
 * Return:
 *  The pointer from malloc()
 */
mcrd_vec* mcrd_alloc_vec(mcrd_int numel);
/*
 * Free all pointers in *x then free x.
 *
 * Args:
 *  x:  a struct to be freed from memory.
 */
void mcrd_free_vec(mcrd_vec* x);
/*
 * Check the sizes of the vector src and dest.
 * If they match, then copy src->c's contents into dest->c.
 * Otherwise throw an error and crash.
 *
 * Args:
 *  x:  a struct to be freed from memory.
 */
void mcrd_copy(mcrd_vec* dest, mcrd_vec* src);

/*
 * Allocate an array of mcrd vectors which are accessed as
 * x[k].c[i] = entry i of vector k. (0 <= k < nvecs)
 * x[k].n    = numel for each k.
 * x[0].c points to the top of a numel by nvecs column major
 * array of all the coordinates. This function only does two mallocs.
 * To free, call
 *   mcrd_free_vec(x)
 *
 * Args:
 *  nvecs:    Number of vectors to be created.
 *  numel:    Number of elemnts for vector.
 *
 * Return:
 *  The pointer from malloc()
 */
mcrd_vec* mcrd_alloc_block(mcrd_int nvecs, mcrd_int numel);

/*
 * Check the sizes of the vectors x, y, z. If they don't match,
 * then throw an error and crash.
 * Otherwise:
 * Linear combination of 2 vectors with n entries.
 * z->c[i] = a*x->c[i] + b*y->c[i]
 * For each i = 0,1,...,n-1
 *
 * Args:
 *  a:  linear combination scalar for x.
 *  x:  vector to be summed
 *  b:  linear combination scalar for y.
 *  y:  vector to be summed
 *  z:  The result of the mathematical sum z = a*x + b*y
 *
 */
void mcrd_axpby(mcrd_flt a, mcrd_vec* x, mcrd_flt b, mcrd_vec* y, mcrd_vec* z);

/*
 * Check the sizes of the vectors x, y, z, w. If they don't match,
 * then throw an error and crash.
 * Otherwise:
 * Linear combination of 3 vectors with n entries.
 * w->c[i] = a*x->c[i] + b*y->c[i] + c*z->c[i]
 * For each i = 0,1,...,n-1
 *
 * Args:
 *  a:  linear combination scalar for x.
 *  x:  vector to be summed
 *  b:  linear combination scalar for y.
 *  y:  vector to be summed
 *  c:  linear combination scalar for z.
 *  z:  vector to be summed
 *  z:  The result of the mathematical sum w = a*x + b*y + c*z
 */
void mcrd_axpbypcz(mcrd_flt a, mcrd_vec* x, mcrd_flt b, mcrd_vec* y, mcrd_flt c,
                   mcrd_vec* z, mcrd_vec* w);

/*
 * Linear combination with variable argument list.
 * variable argument list expected in alternating order, floats then vecs.
 * This summation method uses the compensation summation method described
 * in chapter 4 section 3 of Accuracy and Stability of Numerical Algorithms
 * by Nicholas Higham.
 * Example:
 *      mcrd_lincombo(x,ptrWrk,scalWrk,3, 3.0,z, 4.0,w, -5.0,p);
 * Executes as:
 *          x->c[i] = 3.0*z->c[i] + 4.0*w->c[i] + (-5.0)*p->c[i]
 *
 * Here, we have numTerms = 3, a sum with 3 terms. We also check the sizes
 * of each vector to see if they match x->n. If not, then throw an error
 * and crash.
 * Note: In the example above, the sum is started by setting
 *       x->c[i] = 3.0*z->c[i]
 *       then parsing the remaining arguments in a for loop.
 *       If behavior like:
 *       x = x + y
 *       is desired, then  x must be the first term in the variable argument
 *       list.
 *
 * Args:
 *  x        : The result of a linear combination with numTerms many terms.
 *  ptrWrk   : A mcrd_flt** with at least sizeof(mcrd_flt*)*numTerms memory.
 *             Used for internal pointer arithmetic.
 *  scalWrk  : A mcrd_flt* with sizeof(mcrd_flt)*numTerms memory.
 *             Used for internal pointer arithmetic.
 *  numTerms : Number of terms in linear combination.
 *  ...      : Variable arguments, see above for explanation.
 */
void mcrd_lincombo(mcrd_vec* x, mcrd_flt** ptrWrk, mcrd_flt* scalWrk,
                   mcrd_int numTerms, ...);

/*
 * Check the sizes of the vectors x, y. If they don't match,
 * then throw an error and crash.
 * Otherwise:
 * root mean square difference of x and y. (Weighted 2 norm difference.)
 * Mathematically: If both vectors are of dimension N
 *   sqrt( innerprod(x-y,x-y)/N )
 *
 * Args:
 *  x:  A vector to be diff'd
 *  y:  A vector to be diff'd
 *
 * Return:
 *  Distance beftween x and y  using formula above.
 */
mcrd_flt mcrd_mse(mcrd_vec* x, mcrd_vec* y);

/*
 * root mean square norm of x. (Weighted 2 norm.)
 * Mathematically: If both vectors are of dimension N
 *   sqrt( innerprod(x,y)/N )
 *
 * Args:
 *  x   :  A vector we want to know the norm of.
 *
 * Return:
 *  Distance beftween x and the origin  using formula above.
 */
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
 * x_old    :  mcrd vector of the current step.
 * x_new    :  mcrd vector for order 1 estimate of the next step.
 * vFld_old :  stores the value of vecField(x_old).
 * dt       :  Time step size for ODE solver.
 * vecField :  A function pointer where first argument is input and second
 *             argument is the return value for an ODE's velocity field.
 *             Third variable is the variable argument count, always assumed
 *             zero at runtime. It is assumed that static variables
 *             are used within the implementation of vecField. These are
 *             set using the variable argument count, then not modified
 *             during calls by this routine.
 */
void mcrd_euler_step(mcrd_vec* x_old, mcrd_vec* x_new, mcrd_vec* vFld_old,
                     mcrd_flt dt,
                     void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...));

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
 * x_old    :  mcrd vector of the current step.
 * x_new    :  mcrd vector for order 2 estimate of the next step.
 * vFld_old :  stores the value of vecField(x_old)
 * x_stg    :  stores the value of x_stg
 * dt       :  Time step size for ODE solver.
 * vecField :  A function pointer where first argument is input and second
 *             argument is the return value for an ODE's velocity field.
 *             Third variable is the variable argument count, always assumed
 *             zero at runtime. It is assumed that static variables
 *             are used within the implementation of vecField. These are
 *             set using the variable argument count, then not modified
 *             during calls by this routine.
 *             Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *             variable argument count.
 */
void mcrd_heun_step(mcrd_vec* x_old, mcrd_vec* x_new, mcrd_vec* vFld_old,
                    mcrd_vec* x_stg, mcrd_flt dt,
                    void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...));

/*
 * Use a 2 stage embedded Runge-Kutta method (Heun) to automatically select
 * the a time step satisfying the given absolute tolerance. Used for
 * solving an ODE of the form
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_old       : mcrd vector of the current step.
 * x_new       : mcrd vector estimate of the next step. Guaranteed to be within
 *               abs_tol radius of an order 2 estimate of the time step
 *               pointed to be dt_new_ptr
 * vFld_old    : On input, expected to contain vecField(x_old). The value
 *               is not modified.
 * vFld_new    : On output, contains vecField(x_new). Not expected to have
 *               any particular value on input.
 * x_cmp       : vector containing the comparison step of an order 2 method.
 * dt_old_ptr  : On input, pointer to time step size used to obtain x_old.
 * dt_new_ptr  : On output, pointer to time step size used to obtain x_new.
 * vecField:     A function pointer where first argument is input and second
 *               argument is the return value for an ODE's velocity field.
 *               Third variable is the variable argument count, always assumed
 *               zero at runtime. It is assumed that static variables
 *               are used within the implementation of vecField. These are
 *               set using the variable argument count, then not modified
 *               during calls by this routine.
 *               Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *               variable argument count.
 * err_old_ptr : On input, pointer to the error achieved last time for x_old.
 * err_new_ptr : On output, pointer to the error achieved this time for x_new.
 * nrm_old_ptr : On input, pointer to norm of x_old
 * nrm_new_ptr : On output, pointer to norm of x_new.
 * abs_tol     : The desired local absolute error tolerance.
 * rel_tol     : The desired local relative error tolerance.
 *
 */
void mcrd_o1_autostep(mcrd_vec* x_old, mcrd_vec* x_new, mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new, mcrd_vec* x_cmp, mcrd_flt* dt_old_ptr,
                      mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                      mcrd_flt* err_old_ptr, mcrd_flt* err_new_ptr,
                      mcrd_flt* nrm_old_ptr, mcrd_flt* nrm_new_ptr,
                      mcrd_flt abs_tol, mcrd_flt rel_tol);

/*
 * Use the Heun embedded Runge-Kutta Runge-Kutta method to solve
 * an automous Ordinary Differential Equation of the form.
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_init   : the initial condition to this ODE. This struct is not modified.
 * x_snap   : Pointer to an array of time snapshots of the ODE solution.
 *            These are computed  using a dense output spline formula at the
 *            same order of arracy of this method. The expected memory layout
 *            is specified in mcrd_alloc_block().
 *            To access coorinate j at time k, dereference the x_snap variable
 *            using x_snap[k].c[j]. For each  k, we have x_snap[k].n is
 *            equal to x_init->n.
 * t        : An array of floating point values representing the times at which
 *            we want to know the solution value.
 * t_len    : The number of entries in the array t.
 * vecField : A function pointer where first argument is input and second
 *            argument is the return value for an ODE's velocity field.
 *            Third variable is the variable argument count, always assumed
 *            zero at runtime. It is assumed that static variables
 *            are used within the implementation of vecField. These are
 *            set using the variable argument count, then not modified
 *            during calls by this routine.
 *            Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *            variable argument count.
 * abs_tol  : The desired local absolute error tolerance.
 * rel_tol  : The desired local relative error tolerance.
 * workVec  : A floating point array for storing temporary values. It is
 *            assumed to have 5*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o1(mcrd_vec* x_init, mcrd_vec* x_snap, mcrd_flt* t,
                       mcrd_int t_len,
                       void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                       mcrd_flt abs_tol, mcrd_flt rel_tol, mcrd_flt* workVec);

/*
 * Use the Bogacki-Shampine embedded Runge-Kutta method to automatically select
 * the a time step satisfying the given absolute tolerance. Used for
 * solving an ODE of the form
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_old       : mcrd vector of the current step.
 * x_new       : mcrd vector estimate of the next step. Guaranteed to be within
 *               abs_tol radius of an order 2 estimate of the time step
 *               pointed to be dt_new_ptr
 * vFld_old    : On input, expected to contain vecField(x_old). The value
 *               is not modified.
 * vFld_new    : On output, contains vecField(x_new). Not expected to have
 *               any particular value on input.
 * vFld_stg1   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * vFld_stg2   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * x_cmp       : vector containing the comparison step of an order 2 method.
 * dt_old_ptr  : On input, pointer to time step size used to obtain x_old.
 * dt_new_ptr  : On output, pointer to time step size used to obtain x_new.
 * vecField:     A function pointer where first argument is input and second
 *               argument is the return value for an ODE's velocity field.
 *               Third variable is the variable argument count, always assumed
 *               zero at runtime. It is assumed that static variables
 *               are used within the implementation of vecField. These are
 *               set using the variable argument count, then not modified
 *               during calls by this routine.
 *               Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *               variable argument count.
 * err_old_ptr : On input, pointer to the error achieved last time for x_old.
 * err_new_ptr : On output, pointer to the error achieved this time for x_new.
 * nrm_old_ptr : On input, pointer to norm of x_old
 * nrm_new_ptr : On output, pointer to norm of x_new.
 * abs_tol     : The desired local absolute error tolerance.
 * rel_tol     : The desired local relative error tolerance.
 *
 */
void mcrd_o2_autostep(mcrd_vec* x_old, mcrd_vec* x_new, mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new, mcrd_vec* vFld_stg1,
                      mcrd_vec* vFld_stg2, mcrd_vec* x_cmp,
                      mcrd_flt* dt_old_ptr, mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                      mcrd_flt* err_old_ptr, mcrd_flt* err_new_ptr,
                      mcrd_flt* nrm_old_ptr, mcrd_flt* nrm_new_ptr,
                      mcrd_flt abs_tol, mcrd_flt rel_tol);
/*
 * Use the Bogacki-Shampine embedded Runge-Kutta Runge-Kutta method to solve
 * an automous Ordinary Differential Equation of the form.
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_init   : the initial condition to this ODE. This struct is not modified.
 * x_snap   : Pointer to an array of time snapshots of the ODE solution.
 *            These are computed  using a dense output spline formula at the
 *            same order of arracy of this method. The expected memory layout
 *            is specified in mcrd_alloc_block().
 *            To access coorinate j at time k, dereference the x_snap variable
 *            using x_snap[k].c[j]. For each  k, we have x_snap[k].n is
 *            equal to x_init->n.
 * t        : An array of floating point values representing the times at which
 *            we want to know the solution value.
 * t_len    : The number of entries in the array t.
 * vecField : A function pointer where first argument is input and second
 *            argument is the return value for an ODE's velocity field.
 *            Third variable is the variable argument count, always assumed
 *            zero at runtime. It is assumed that static variables
 *            are used within the implementation of vecField. These are
 *            set using the variable argument count, then not modified
 *            during calls by this routine.
 *            Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *            variable argument count.
 * abs_tol  : The desired local absolute error tolerance.
 * rel_tol  : The desired local relative error tolerance.
 * workVec  : A floating point array for storing temporary values. It is
 *            assumed to have 7*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o2(mcrd_vec* x_init, mcrd_vec* x_snap, mcrd_flt* t,
                       mcrd_int t_len,
                       void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                       mcrd_flt abs_tol, mcrd_flt rel_tol, mcrd_flt* workVec);

/*
 * Use the DOPRI5 Runge-Kutta method to automatically select
 * the a time step satisfying the given absolute tolerance. Used for
 * solving an ODE of the form
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_old       : mcrd vector of the current step.
 * x_new       : mcrd vector estimate of the next step. Guaranteed to be within
 *               abs_tol radius of an order 2 estimate of the time step
 *               pointed to be dt_new_ptr
 * vFld_old    : On input, expected to contain vecField(x_old). The value
 *               is not modified.
 * vFld_new    : On output, contains vecField(x_new). Not expected to have
 *               any particular value on input.
 * vFld_stg1   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * vFld_stg2   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * vFld_stg3   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * vFld_stg4   : A work variable used for computing Runge-Kutta stages.
 *               Expected to have the same size as vFld_old and vFld_new
 * x_cmp       : vector containing the comparison step of an order 4 method.
 * dt_old_ptr  : On input, pointer to time step size used to obtain x_old.
 * dt_new_ptr  : On output, pointer to time step size used to obtain x_new.
 * vecField:     A function pointer where first argument is input and second
 *               argument is the return value for an ODE's velocity field.
 *               Third variable is the variable argument count, always assumed
 *               zero at runtime. It is assumed that static variables
 *               are used within the implementation of vecField. These are
 *               set using the variable argument count, then not modified
 *               during calls by this routine.
 *               Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *               variable argument count.
 * err_old_ptr : On input, pointer to the error achieved last time for x_old.
 * err_new_ptr : On output, pointer to the error achieved this time for x_new.
 * nrm_old_ptr : On input, pointer to norm of x_old
 * nrm_new_ptr : On output, pointer to norm of x_new.
 * abs_tol     : The desired local absolute error tolerance.
 * rel_tol     : The desired local relative error tolerance.
 *
 */
void mcrd_o4_autostep(mcrd_vec* x_old, mcrd_vec* x_new, mcrd_vec* vFld_old,
                      mcrd_vec* vFld_new, mcrd_vec* vFld_stg1,
                      mcrd_vec* vFld_stg2, mcrd_vec* vFld_stg3,
                      mcrd_vec* vFld_stg4, mcrd_vec* x_cmp,
                      mcrd_flt* dt_old_ptr, mcrd_flt* dt_new_ptr,
                      void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                      mcrd_flt* err_old_ptr, mcrd_flt* err_new_ptr,
                      mcrd_flt* nrm_old_ptr, mcrd_flt* nrm_new_ptr,
                      mcrd_flt abs_tol, mcrd_flt rel_tol);

/*
 * Use the DOPRI5 embedded Runge-Kutta Runge-Kutta method to solve
 * an automous Ordinary Differential Equation of the form.
 *
 *  dx
 * ----  = vecField(x)
 *  dt
 *
 * x_init   : the initial condition to this ODE. This struct is not modified.
 * x_snap   : Pointer to an array of time snapshots of the ODE solution.
 *            These are computed  using a dense output spline formula at the
 *            same order of arracy of this method. The expected memory layout
 *            is specified in mcrd_alloc_block().
 *            To access coorinate j at time k, dereference the x_snap variable
 *            using x_snap[k].c[j]. For each  k, we have x_snap[k].n is
 *            equal to x_init->n.
 * t        : An array of floating point values representing the times at which
 *            we want to know the solution value.
 * t_len    : The number of entries in the array t.
 * vecField : A function pointer where first argument is input and second
 *            argument is the return value for an ODE's velocity field.
 *            Third variable is the variable argument count, always assumed
 *            zero at runtime. It is assumed that static variables
 *            are used within the implementation of vecField. These are
 *            set using the variable argument count, then not modified
 *            during calls by this routine.
 *            Evaluated as vecField(x_old,x_new,0). 3rd argument is the
 *            variable argument count.
 * abs_tol  : The desired local absolute error tolerance.
 * rel_tol  : The desired local relative error tolerance.
 * workVec  : A floating point array for storing temporary values. It is
 *            assumed to have 9*x_init->n many mcrd_flts allocated.
 */
void mcrd_ode_solve_o4(mcrd_vec* x_init, mcrd_vec* x_snap, mcrd_flt* t,
                       mcrd_int t_len,
                       void (*vecField)(mcrd_vec*, mcrd_vec*, int, ...),
                       mcrd_flt abs_tol, mcrd_flt rel_tol, mcrd_flt* workVec);

#endif
