# virasoro
Numerical computation of degenerate Virasoro block coefficients in CFT_2.

To compile, type "make" into a terminal in this directory. You will need GCC (or another compiler if you change the makefile) and the GMP arbitrary-precision math library, which you can download from gmplib.org.

Some interior parameters are included in static variables initialized at the top of virasoro.h: the number of threads used in the slow part of the computation (default 8, appropriate for 4-core processors) and the floating-point precision of the numbers used (default 512, which allows low-error calculations at least up to order 1000). It should still compile and run fine with any reasonable values of these numbers.
