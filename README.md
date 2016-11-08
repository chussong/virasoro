# virasoro
Numerical computation of Virasoro block coefficients in CFT_2.

To compile, type "make" into a terminal in this directory. You will need GCC (or another compiler if you change the makefile) and the GMP arbitrary-precision math library, which you can download from gmplib.org. GMP must be built with C++ support by using "./configure --enable-cxx" before invoking its make; prebuilt versions of GMP are available for some systems. If you get an error saying libgmp.so.4 is missing, it means you don't have GMP installed correctly; if you get one saying libgmpxx.so.4 is missing, it means you have installed GMP correctly but without C++ support.

There are two ways to invoke runs. For a single run, simply pass the parameters c, hl, hh, hp, and qmax. For example, to get a run with c=30, hl=1, hh=3, hp=0, qmax=1000, you would type:  
./virasoro 30 1 3 0 1000  
To do multiple runs at once, you can input a "runfile". Instead of five parameters, you give one: the name of the file. For example:  
./virasoro 1to30.txt  
will tell the program to read 1to30.txt and perform all the runs shown in it. A runfile should be plain text, with one run per line. For example, a runfile containing the following two lines:  
30 1 3 0 1000  
35 1 3 0 1000  
would run the program twice, once with c=30 and once with c=35. Parameters can be separated by single spaces, commas, or semicolons. Multiple runs can be described in a batch using Mathematica-like syntax: for instance, the following line:  
{26,34,2} 1 3 {0,5,1} 1000  
would expand to running with hp=0,1,2,3,4,5 for each of c=26,28,30,32,34, a total of 30 order-1000 runs. This would most likely finish in under an hour. The program automatically combines runs which differ only in hp and/or qmax in order to save time.

There are two known seemingly-plausible input categories which will cause the program to fail. First, c between 1 and 25 will result in a complex b^2, which it currently can not handle. Second, if b^2 or 1/b^2 is an integer, then divisions by 0 appear; this occurs for example when c=25, 28, 33, 38.5, 44.2, 50, etc., and runs with these values of c are skipped automatically up to b^2 = 20.

There are two options available, which can be placed either before or after the arguments. They are:  
-m		|		give output as a single Mathematica vector suitable for return from RunThrough  
-c		|		print output to the console but do not write a file

To call this function from Mathematica, do the following:  
SetDirectory["\<directory where virasoro is\>"];  
vec = RunThrough["./virasoro -m \<c\> \<hl\> \<hh\> \<hp\> \<maxOrder\>",""]; or "./virasoro -m \<name of runfile\>"  
Either one of these calls will fill vec with a Mathematica vector containing parameters used (in odd entries) and coefficients (in even entries).

Some interior parameters are included in static variables initialized at the top of virasoro.cpp: the number of threads used in the slow part of the computation (default 8, appropriate for 4-core processors) and the floating-point precision of the numbers used (default 512, which allows low-error calculations at least up to order 1000). It should still compile and run fine with any reasonable values of these numbers.

If you want to read the source code, it's almost certainly better to start with the deprecated 128-bit version in the backup folder. The algorithm is the same, and the __float128s use standard mathematical syntax, unlike the GMP numbers in the current version.
