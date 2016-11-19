#ifndef RUNFILE_H_ 
#define RUNFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <gmpxx.h>
#include "virasoro.h"

class Runfile_c{
	const std::string filename;
	std::vector<std::string> lines;

	public:
		std::vector<std::vector<mpf_class>> runs;
		std::vector<int> maxOrders;

		Runfile_c(const char* filename);
		Runfile_c(const std::string filename);
		
		int ReadRunfile();

		int Expand();

		int ExpandBraces();

		std::tuple<mpf_class, mpf_class, mpf_class> ParseBraces(std::string firstHalf, std::string insideBraces);

		int ExpandRelativeEqns();
			
		std::tuple<mpf_class, int> ParseRelativeEqn(std::string equation, std::string relTo);

		mpf_class RelativeMPF(std::string firstHalf, std::string equation);

		int ReadMPF(mpf_class& output, FILE* runfile);

		int ReadMaxOrder(FILE* runfile);

		int RunCompare(std::vector<mpf_class> run1, std::vector<mpf_class> run2);
};

#endif
