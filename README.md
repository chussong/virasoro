# virasoro
A program for numerical computation of Virasoro block coefficients in CFT<sub>2</sub>.  

Copyright 2016-2017 Charles Hussong  

Project homepage:	https://github.com/chussong/virasoro  
Contact email:		chussong@jhu.edu  

## License

Virasoro is licensed under the GNU General Public License, enclosed as LICENSE.txt. It is free software, and you can redistribute and/or modify it under the terms of the GPL, either version 3 of the License or any later version.  

This work makes extensive use of [GNU MP](https://gmplib.org/), its floating-point extension [GNU MPFR](http://www.mpfr.org/), and the complex-number extension of MPFR, [GNU MPC](http://www.multiprecision.org/index.php?prog=mpc). The new container for MPC, mpcomplex.h, is based on [mpreal](http://www.holoborodko.com/pavel/mpfr/), a container for MPFR.  

## Installation

To compile, get a terminal in this directory and type ./configure, then make, then sudo make install. See the file called INSTALL for build options. You will need working C and C++ compilers as well as the following libraries installed in this order:  
0. Mac users will need to run xcode-select --install to get the Command Line Tools.  
1. GMP arbitrary-precision math library from [here](https://gmplib.org).  
2. GNU MPFR from [here](http://www.mpfr.org).  
3. GNU MPC from [here](http://www.multiprecision.org/index.php?prog=mpc).  
4. Linux users who want to build the Mathematica plugin will need to install
libuuid, available on apt as uuid-dev and yum as libuuid-devel. If you have 
Mathematica, the installer will automatically try to build the plugin, but you
can manually override this with ./configure --disable-wstp.  
Prebuilt versions of GMP, MPFR, and MPF are available for some systems (they're on the Ubuntu repositories, anyway), but in my brief testing the prebuilt ones were about half as fast.  

## Usage

You can run this program either from a terminal or Mathematica, and can either do runs one at a time or in a large batch. For a single run, simply pass the parameters c, hl, hh, h, and qmax. For example, to get a run with c=30, hl=1, hh=3, h=0, qmax=1000, you would type:  
virasoro 30 1 3 0 1000  

To do multiple runs at once, you can input a "runfile". Instead of five parameters, you give one: the name of the file. For example:  
virasoro testruns.txt  
will tell the program to read testruns.txt and perform all the runs shown in it. A runfile should be plain text, with one run per line. For example, a runfile containing the following two lines:  
30 1 3 0 1000  
35 1 3 0 1000  
would run the program twice, once with c=30 and once with c=35. Parameters can be separated by single spaces, commas, or semicolons. A larger example runfile called testruns.txt is included in the repository.  

Multiple runs can be described in a batch on one line using Mathematica-like syntax: for instance, the following line:  
{26,34,2} 1 3 {0,5,1} 1000  
would expand to running with h=0,1,2,3,4,5 for each of c=26,28,30,32,34, a total of 30 order-1000 runs. This would most likely finish in under an hour. The program automatically combines runs which differ only in h and/or qmax in order to save time.  
It's also possible to specify parameters in terms of others, which is quite useful for batch runs. For example, if you wanted to do a lot of runs with hl = c/25 and hh = c/11, you could just write the following:  
{26,34,0.2} c/25 c/11 0 1000  
This is still not very robust, so you have to be careful: the parameters must be immediately followed or preceded by an explicit operation, so c/25 and 25\*c are legal but c / 25 and 25c are not. Parameters can also only be defined relative to those which appear to their left, so hh can be set to 2.3\*hl but hl can not be set to hh/2.3. +, -, \*, and / are supported, as well as a parameter just sitting by itself, e.g. hl.  
Complex parameters can be entered as (realPart imPart), e.g. 1 + i is (1 1). **If you want to invoke runs using parentheses or brackets from the command line, they must be included in quotation marks, like so:**  
virasoro "{10,20,2}" "(1 1)" 2\*hl 0 1000  

There are a couple of parameter restrictions to be aware of. First, if c is between 1 and 25, the computation will have to be done using complex numbers, which is about 8 times slower. Second, if b^2 or 1/b^2 is a rational number with relatively small numerator or denominator, it will cause the internal Hmn and/or Amn to diverge. If such a c is detected, the maxOrder will automatically be adjusted down to the highest safe value.  

## Options

There are six options available, which can be placed either before or after the arguments. They are:  

| Option | Description |
| ------ | ----------- |
| -m | give output as a single \[M\]athematica vector suitable for return from RunThrough |
| -c | print output to the \[c\]onsole but do not write a file |
| -p# | set numerical \[p\]recision to # (any positive integer, default would be -p768) |
| -t# | set number of \[t\]hreads to # (any positive integer, default would be -t8) |
| -b | first provided value is \[b\] instead of c |
| -bb | first provided value is \[b^2\] instead of c |

There is also a configuration file, created at ~/.config/virasoro\_defaults.txt, which contains the default values of these options so you can change them there as well. The command line options override the entries in virasoro\_defaults.txt. The entries mean the following:

| Option | Description |
| ------ | ----------- |
| maxThreads | Maximum number of threads to use in the slow part of the calculation. Default 8, appropriate for a 4-core system. |
| precision | Number of bits of precision used to represent all numbers internally. Default is 768, which is good enough to get most parameter ranges to order 1000. |
| tolerance | This is the level below which numbers will be treated as zero when writing output (full precision is always kept internally). Default is 10^-100. |
| showProgressBar | Shows a progress bar for the slow part of the calculation. Automatically suppressed if -m is given in run commands. This defaults to true; set it to false if you don't like it or really need those precious milliseconds. |

## Mathematica

If Mathematica is detected, a Mathematica package will also be installed. On Linux this is saved to ~/.Mathematica/Applications/Virasoro.m; on Mac OS it's saved to ~/Library/Mathematica/Applications/Virasoro.m.  

To run this program from Mathematica, load Virasoro.m by calling Needs["Virasoro`"] \(it should be located automatically\) and call VRun, like one of the following:  
results = VRun[30, 1, 3, 0, 1000];  
results = VRun[runfile.txt];  
If VWSTP was built and installed correctly, this will fill "results" with a Mathematica vector of vectors containing parameters used (in odd entries) and coefficients (in even entries). This can then be plotted with VPlot[results]; for explanations of the many options for VPlot, call ?VPlot.  

You can also read in the results of runs invoked from the command line using VRead["filename"], which will give the same output as if they were invoked with VRun. If you do runs within Mathematica and want to keep the results, call VWrite[results, "filename"] to save them as filename. Other useful functions in Virasoro.m are VPlotCoeffs, VConvByOrder, and VConvByTL; once you've run Needs["Virasoro`"], descriptions of these functions are all available from ?FunctionName.  

If you would like to automatically load Virasoro when Mathematica starts, add a line containing Needs["Virasoro`"] to ~/.Mathematica/Kernel/init.m. Warning: running this takes a few seconds, so you probably don't want to do it. To stop the automatic loading, simply delete the line from init.m.  
