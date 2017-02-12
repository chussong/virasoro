#ifndef RUNFILE_H_ 
#define RUNFILE_H_

#include <iostream>
#include <fstream>
#include <algorithm>	// std::swap
#include <string>
#include <vector>
#include <tuple>
#ifdef HAVE_WSTP
#include "wstp.h"
#endif
#include "mpreal.h"
#include "mpcomplex.h"
//#include "config.h" !! is this actually necessary?

namespace virasoro {

class Runfile_c{
	public:
		std::string filename;
		std::vector<std::string> lines;

		std::vector<std::vector<std::complex<mpfr::mpreal>>> runs;
		std::vector<int> maxOrders;

		Runfile_c(const char* filename);
		Runfile_c(const std::string filename);
		Runfile_c(const std::vector<std::string> line);

		int NumberOfRuns();

		int ReadRunfile();

		int Expand();
		int ExpandBraces(const int param);
		std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> ParseBraces(std::string insideBraces);
		std::tuple<size_t, size_t> FindNthParameter(const std::string line, const int param);

		int ExpandRelativeEqns(const int param);
		std::tuple<std::complex<mpfr::mpreal>, int> ParseRelativeEqn(std::string equation, std::string relTo);
		std::complex<mpfr::mpreal> RelativeMPF(std::string firstHalf, std::string equation);
		std::string FindBaseNumber(std::string sourceString, const int paramNumber);

		int RunCompare(std::vector<std::complex<mpfr::mpreal>> run1, std::vector<std::complex<mpfr::mpreal>> run2);
};

} // namespace virasoro
#endif
