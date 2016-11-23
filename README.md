# virasoro
Numerical computation of Virasoro block coefficients in CFT_2.

To compile, type "make" into a terminal in this directory. You will need GCC (or another compiler if you change the makefile) and the following dependencies:
GMP arbitrary-precision math library from [here](gmplib.org).  
GNU MPFR from [here](mpfr.org).  
GNU MPC from [here](multiprecision.org/index.php?prog=mpc&page=home).  
Prebuilt versions of these are available for some systems (they're on the Ubuntu repositories, anyway); if you build it yourself, GMP must be compiled with C++ support by using "./configure --enable-cxx" before invoking its make. If you get an error saying libgmp.so.4 is missing, it means you don't have GMP installed correctly; if you get one saying libgmpxx.so.4 is missing, it means you have installed GMP correctly but without C++ support.  

There are two ways to invoke runs. For a single run, simply pass the parameters c, hl, hh, hp, and qmax. For example, to get a run with c=30, hl=1, hh=3, hp=0, qmax=1000, you would type:  
./virasoro 30 1 3 0 1000  

To do multiple runs at once, you can input a "runfile". Instead of five parameters, you give one: the name of the file. For example:  
./virasoro testruns.txt  
will tell the program to read 1to30.txt and perform all the runs shown in it. A runfile should be plain text, with one run per line. For example, a runfile containing the following two lines:  
30 1 3 0 1000  
35 1 3 0 1000  
would run the program twice, once with c=30 and once with c=35. Parameters can be separated by single spaces, commas, or semicolons.  

Multiple runs can be described in a batch using Mathematica-like syntax: for instance, the following line:  
{26,34,2} 1 3 {0,5,1} 1000  
would expand to running with hp=0,1,2,3,4,5 for each of c=26,28,30,32,34, a total of 30 order-1000 runs. This would most likely finish in under an hour. The program automatically combines runs which differ only in hp and/or qmax in order to save time.  
It's also possible to specify parameters in terms of others, which is quite useful for batch runs. For example, if you wanted to do a lot of runs with hl = c/25 and hh = c/11, you could just write the following:  
{26, 34, 0.2} c/25 c/11 0 1000  

This is still not very robust, so you have to be careful: the parameters must be immediately followed or preceded by an explicit operation, so c/25 and 25\*c are legal but c / 25 and 25c are not. Parameters can also only be defined relative to those which appear to their left, so hh can be set to 2.3\*hl but hl can not be set to hh/2.3. +, -, \*, and / are supported. **The bracket and arithmetic syntax only works within runfiles, not from the terminal.**  

There are a couple of parameter restrictions to be aware of. First, if c is between 1 and 25, the computation will have to be done using complex numbers, which is about 8 times slower; the output is also not ideally set up yet. Second, if b^2 or 1/b^2 is a rational number with relatively small numerator or denominator, it will cause the internal Hmn and/or Amn to diverge. If such a c is detected, the maxOrder will automatically be adjusted down to the highest safe value.  

There are six options available, which can be placed either before or after the arguments. They are:  

| Option | Description |
| ------ | ----------- |
| -m | give output as a single \[M\]athematica vector suitable for return from RunThrough |
| -c | print output to the \[c\]onsole but do not write a file |
| -p# | set numerical \[p\]recision to # (any positive integer, default would be -p512) |
| -t# | set number of \[t\]hreads to # (any positive integer, default would be -t8) |
| -b | first provided value is \[b\] instead of c |
| -bb | first provided value is \[b^2\] instead of c |

To call this function from Mathematica, do the following:  
SetDirectory["\<directory where virasoro is\>"];  
vec = RunThrough["./virasoro -m \<c\> \<hl\> \<hh\> \<hp\> \<maxOrder\>",""]; or "./virasoro -m \<name of runfile\>"  
Either one of these calls will fill vec with a Mathematica vector containing parameters used (in odd entries) and coefficients (in even entries).  

Some interior parameters are included in static variables initialized at the top of virasoro.cpp: the number of threads used in the slow part of the computation (default 8, appropriate for 4-core processors) and the floating-point precision of the numbers used (default 512, which allows low-error calculations at least up to order 1000). It should still compile and run fine with any reasonable values of these numbers.
