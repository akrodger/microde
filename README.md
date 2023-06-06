# microde by Bram Rodgers
An extremely bare-bones ODE toolkit written in C with MIT license. Please
share this code everywhere you can so nobody needs to implement this again.
The entire library is less than 1400 lines when you include the documentation.

# What is this?
microde (Micro ODE) is an implementation of several Runge-Kutta solvers from
order 1 to 5. Automatic time step selection is implemented for the schemes:

- Euler-Heun, an embedded order (1,2) method. (Order 1 due to linear interp.)
- Bogacki-Shampine, an embedded (2,3) method. (Order 2 due to quad interp.)
- Dormand-Prince, an embedded (5,4) method. (Order 4 due to quartic interp.)

(See Solving Ordinary Differential Equations I: Nonstiff Problems
by Hairer, Norsett, and Wanner.)

The code is written to be entirely stand-alone. The only requirement outside
of the C standard library is that the vectors be contiguously allocated arrays.
The integer type and floating point type are customizable via the microde.h
header file. OpenMP is available and may be enabled for all addition operations
through the use of a #define guard. The minimum  number of vector field
operations are done at every time step, making this code have near the minimum
operation count. Compactness and portability are favored wherever possible. The
intent is to make it just as easily compiled on microcontrollers with FPUs
as with shared memory supercomputers.

# How do I use this?
This is intended to be a small, portable, back-end library you can place into
your source files when you need to solve some ODEs as a component of your code.
Documentation for every function is clearly defined in src/microde.h.
The core features of automatic time stepping are in the funtions:

mcrd_ode_solve_o1(), mcrd_ode_solve_o2(), mcrd_ode_solve_o4()

which are the order 1, order 2, and order 4 solvers respectively.

# How do I compile this?
The same way any C header + source pair are compiled. Make sure that
C standard library math is linked as well.
If OpenMP prarllelism is desired, compile with the extra flag:
"-D SHMEM_PARA_MICRODE"
and add the proper OpenMP compiler options at object building, linking, etc.

# I have a problem with (x,y,z)!
File an issue, though I don't guarantee that I will have time to address it.
If it is a big problem and you want to add some massive feature I don't have,
then just fork my code.
