/*
 * Microde DOPRI(5,4) with auto step selection test file by Bram Rodgers.
 * Original Draft Dated: 15, July 2022
 */

/*
 * Macros and Includes go here: (Some common ones listed)
 */
#include <stdio.h>
#include <stdlib.h>

#include "../src/microde.h"
#if defined(SHMEM_PARA_MICRODE)
#include "omp.h"
#endif
/*
 * Function declarations go here:
 */

void vecField(mcrd_vec* x, mcrd_vec* dxdt, int argc, ...);

mcrd_flt* linspace(mcrd_flt left, mcrd_flt right, mcrd_int n);
/*
 * Template main:
 */
int main(int argc, char** argv) {
  char* pEnd = NULL;
  mcrd_flt pi = acos(-1.0);
  mcrd_flt t_final = 1.0;
  mcrd_flt abs_tol = 5e-14;
  mcrd_flt rel_tol = 5e-14;
  mcrd_int Nt;
  mcrd_int numel = 100;
  if (argc > 2) {
    Nt = (mcrd_int)strtol(argv[1], &pEnd, 0);
    numel = (mcrd_int)strtol(argv[2], &pEnd, 0);
  } else {
    Nt = 100;
  }
  mcrd_int k = 0;
  mcrd_flt a = 1.0;
  mcrd_flt dt = t_final;
  mcrd_flt* t = linspace(0.0, t_final, Nt + 1);
  mcrd_vec* x_init = mcrd_alloc_vec(numel);
  mcrd_flt* work = (mcrd_flt*)malloc(sizeof(mcrd_flt) * 9 * numel);
  mcrd_vec* x_snap = mcrd_alloc_block(Nt + 1, numel);
  mcrd_vec* sten = mcrd_alloc_vec(5);
  mcrd_flt gc1[2];
  mcrd_flt gc2[2];
  dt = dt / Nt;
#if defined(SHMEM_PARA_MICRODE)
  double wtime;
#endif
  sten->c[0] = (numel)*1.0 / (12.0);
  sten->c[1] = -(numel + 0) * 2.0 / (3.0);
  sten->c[2] = 0.0;
  sten->c[3] = (numel + 0) * 2.0 / (3.0);
  sten->c[4] = -(numel)*1.0 / (12.0);
  vecField(NULL, NULL, 6, a, 0, sten->c, sten->n, &(gc1[0]), &(gc2[0]));
  for (k = 0; k < numel; k++) {
    x_init->c[k] = exp(sin((2 * pi * k) / numel));
  }
#if defined(SHMEM_PARA_MICRODE)
  wtime = omp_get_wtime();
#endif
  mcrd_ode_solve_o4(x_init, x_snap, t, Nt + 1, &vecField, abs_tol, rel_tol,
                    work);
#if defined(SHMEM_PARA_MICRODE)
  wtime = omp_get_wtime() - wtime;
  printf("\nComp Time = %le", wtime);
#endif
  printf("\nMSE       = %2.16le", (double)mcrd_mse(x_init, &(x_snap[Nt])));
  printf("\n");
  vecField(NULL, NULL, -1);
  mcrd_free_vec(x_init);
  mcrd_free_vec(sten);
  free(work);
  free(t);
  mcrd_free_vec(x_snap);
  return 0;
}

void vecField(mcrd_vec* x, mcrd_vec* dxdt, int argc, ...) {
  static mcrd_flt a;
  static mcrd_int feval;
  static mcrd_flt* sten;
  static mcrd_int slen;
  static mcrd_flt* gc1;
  static mcrd_flt* gc2;
  mcrd_int sm1o2, i, j, slen_bot, N, loopStart, loopEnd, thNum, thTot,
      chunkSize, chunkRem;
  thNum = 0;
  thTot = 1;
  va_list arglist;
  if (argc == -1) {
    printf("Num function evals = %ld\n", (long)feval);
  } else if (argc == 1) {
    va_start(arglist, argc);
    a = va_arg(arglist, mcrd_flt);
  } else if (argc == 6) {
    va_start(arglist, argc);
    a = va_arg(arglist, mcrd_flt);
    feval = va_arg(arglist, mcrd_int);
    sten = va_arg(arglist, mcrd_flt*);
    slen = va_arg(arglist, mcrd_int);
    gc1 = va_arg(arglist, mcrd_flt*);
    gc2 = va_arg(arglist, mcrd_flt*);
    if (slen % 2 == 0) {
      slen = slen - 1;
    }
  } else if (argc == 0) {
    feval = feval + 1;
    sm1o2 = (slen - 1) / 2;
    for (i = 0; i < sm1o2; i++) {
      gc1[i] = x->c[x->n - sm1o2 + i];
      gc2[i] = x->c[i];
    }
#if defined(SHMEM_PARA_MICRODE)
#pragma omp parallel default(shared) private( \
    i, j, slen_bot, chunkSize, N, loopStart, loopEnd, chunkRem, thNum, thTot)
    {
      thTot = omp_get_num_threads();
      thNum = omp_get_thread_num();
#endif
      N = x->n;
      chunkSize = N / thTot;
      chunkRem = N % thTot;
      if (thNum < chunkRem) {
        loopStart = (chunkSize + 1) * thNum;
        loopEnd = loopStart + chunkSize + 1;
      } else {
        loopStart = (chunkSize + 1) * chunkRem + chunkSize * (thNum - chunkRem);
        loopEnd = loopStart + chunkSize;
      }
      for (i = loopStart; i < loopEnd; i++) {
        dxdt->c[i] = 0.0;
        if (i < sm1o2) {
          for (j = 0; j < (sm1o2)-i; j++) {
            dxdt->c[i] += a * sten[j] * gc1[j + i];
          }
          for (j = 0; j <= i + (sm1o2); j++) {
            dxdt->c[i] += a * sten[j + (sm1o2)-i] * x->c[j];
          }
        } else if (i < x->n - (sm1o2)) {
          for (j = 0; j < slen; j++) {
            dxdt->c[i] += a * sten[j] * x->c[i + j - sm1o2];
          }
        } else {  // bottom layer of gohst cells
          for (j = i - (sm1o2); j < N; j++) {
            dxdt->c[i] += a * sten[j + (sm1o2)-i] * x->c[j];
          }
          slen_bot = sm1o2 - (N - i) + 1;
          for (j = 0; j < slen_bot; j++) {
            dxdt->c[i] += a * sten[slen - slen_bot + j] * gc2[j];
          }
        }
      }
#if defined(SHMEM_PARA_MICRODE)
    }
#endif
  }
}

mcrd_flt* linspace(mcrd_flt left, mcrd_flt right, mcrd_int n) {
  mcrd_int k;
  mcrd_flt* x = (mcrd_flt*)malloc(sizeof(mcrd_flt) * n);
  for (k = 0; k < n; k++) {
    x[k] = left + ((k * (right - left)) / (n - 1));
  }
  return x;
}

/*
//Coefs for an order 6 stencil.
sten->c[0] = -(numel+0)*1.0/(60.0);
sten->c[1] =  (numel+0)*3.0/(20.0);
sten->c[2] = -(numel+0)*3.0/(4.0);
sten->c[3] =  0.0;
sten->c[4] =  (numel+0)*3.0/(4.0);
sten->c[5] = -(numel+0)*3.0/(20.0);
sten->c[6] =  (numel+0)*1.0/(60.0);
*/
