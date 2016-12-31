# virasoro
Numerical computation of Virasoro block coefficients in CFT<sub>2</sub>.  

To compile, type "make" into a terminal in this directory. You will need GCC (or another compiler if you change the makefile) and the following dependencies:  
GMP arbitrary-precision math library from [here](https://gmplib.org).  
GNU MPFR from [here](http://www.mpfr.org).  
GNU MPC from [here](http://www.multiprecision.org/index.php?prog=mpc).  
Prebuilt versions of these are available for some systems (they're on the Ubuntu repositories, anyway), but in my brief testing the prebuilt ones were about half as fast. If you build everything yourself, GMP must be compiled with C++ support by using "./configure --enable-cxx" before invoking its make. If you get an error saying libgmp.so.4 is missing, it means you don't have GMP installed correctly; if you get one saying libgmpxx.so.4 is missing, it means you have installed GMP correctly but without C++ support.  
Linux users who want to build the Mathematica plugin will likely also need to install libuuid, available on apt as uuid-dev and yum as libuuid-devel.  

There are two ways to invoke runs. For a single run, simply pass the parameters c, hl, hh, hp, and qmax. For example, to get a run with c=30, hl=1, hh=3, hp=0, qmax=1000, you would type:  
./virasoro 30 1 3 0 1000  

To do multiple runs at once, you can input a "runfile". Instead of five parameters, you give one: the name of the file. For example:  
./virasoro testruns.txt  
will tell the program to read testruns.txt and perform all the runs shown in it. A runfile should be plain text, with one run per line. For example, a runfile containing the following two lines:  
30 1 3 0 1000  
35 1 3 0 1000  
would run the program twice, once with c=30 and once with c=35. Parameters can be separated by single spaces, commas, or semicolons. A larger example runfile called testruns.txt is included in the repository.  

Multiple runs can be described in a batch on one line using Mathematica-like syntax: for instance, the following line:  
{26,34,2} 1 3 {0,5,1} 1000  
would expand to running with hp=0,1,2,3,4,5 for each of c=26,28,30,32,34, a total of 30 order-1000 runs. This would most likely finish in under an hour. The program automatically combines runs which differ only in hp and/or qmax in order to save time.  
It's also possible to specify parameters in terms of others, which is quite useful for batch runs. For example, if you wanted to do a lot of runs with hl = c/25 and hh = c/11, you could just write the following:  
{26,34,0.2} c/25 c/11 0 1000  
This is still not very robust, so you have to be careful: the parameters must be immediately followed or preceded by an explicit operation, so c/25 and 25\*c are legal but c / 25 and 25c are not. Parameters can also only be defined relative to those which appear to their left, so hh can be set to 2.3\*hl but hl can not be set to hh/2.3. +, -, \*, and / are supported.  
Complex parameters can be entered as (realPart imPart), e.g. 1 + i is (1 1). **If you want to invoke runs using parentheses or brackets from the command line, they must be included in quotation marks, like so:**  
./virasoro "{10,20,2}" "(1 1)" 2\*hl 0 1000  

There are a couple of parameter restrictions to be aware of. First, if c is between 1 and 25, the computation will have to be done using complex numbers, which is about 8 times slower. Second, if b^2 or 1/b^2 is a rational number with relatively small numerator or denominator, it will cause the internal Hmn and/or Amn to diverge. If such a c is detected, the maxOrder will automatically be adjusted down to the highest safe value.  

There are six options available, which can be placed either before or after the arguments. They are:  

| Option | Description |
| ------ | ----------- |
| -m | give output as a single \[M\]athematica vector suitable for return from RunThrough |
| -c | print output to the \[c\]onsole but do not write a file |
| -p# | set numerical \[p\]recision to # (any positive integer, default would be -p512) |
| -t# | set number of \[t\]hreads to # (any positive integer, default would be -t8) |
| -b | first provided value is \[b\] instead of c |
| -bb | first provided value is \[b^2\] instead of c |

There is also a configuration file, config.txt, which contains the default values of these options so you can change them there as well. The command line options override the entries in config.txt. The entries mean the following:

| Option | Description |
| ------ | ----------- |
| maxThreads | Maximum number of threads to use in the slow part of the calculation. Default 8, appropriate for a 4-core system. |
| precision | Number of bits of precision used to represent all numbers interally. Default is 512, which is good enough to get most parameter ranges to order 1000. |
| tolerance | This is the level below which numbers will be treated as zero when checking, for example, if a complex number is purely real. Default 10^-20, though that's probably aggressive. |
| showProgressBar | Shows a progress bar for the slow part of the calculation. Automatically suppressed if -m is given in run commands. Default true, set to false if you don't like it or really need those precious milliseconds. |

To call this function from Mathematica, load the library Virasoro.m like in the provided vwstp\_test.nb and call VRun, like one of the following:  
results = VRun[30, 1, 3, 0, 1000];  
results = VRun[runfile.txt];  
This will fill results with a Mathematica vector of vectors containing parameters used (in odd entries) and coefficients (in even entries). This can then be plotted with VPlot[results, 1, 1, 0.99]; for explanations of the parameters to VPlot, call ?VPlot.
