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
		std::vector<std::vector<std::string>> entries;

		std::vector<std::vector<std::complex<mpfr::mpreal>>> runs;
		std::vector<int> maxOrders;

		Runfile_c(const char* filename);
		Runfile_c(const std::string filename);
		Runfile_c(const std::vector<std::string> line);

		int NumberOfRuns();

		int ReadRunfile();

		int Expand();
		int ExpandBraces(const int param);
		void ExpandRelativeEqns(const int param);
		std::complex<mpfr::mpreal> RelativeMPF(int lineNum, std::string equation);

		static std::vector<std::vector<std::string>> ParseLines(const std::vector<std::string>& lines);
		static std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> ParseBraces(std::string insideBraces);
		static std::tuple<std::complex<mpfr::mpreal>, int> ParseRelativeEqn(std::string equation, std::string relTo);
		static std::tuple<std::vector<std::vector<std::complex<mpfr::mpreal>>>,std::vector<int>> NumericizeRuns(const std::vector<std::vector<std::string>>& entries);
		static int RunCompare(const std::vector<std::complex<mpfr::mpreal>>& run1, const std::vector<std::complex<mpfr::mpreal>>& run2);
		static void CrunchDuplicateRuns(std::vector<std::vector<std::complex<mpfr::mpreal>>>& runs, std::vector<int>& maxOrders);
};

} // namespace virasoro
#endif
