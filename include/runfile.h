#ifndef RUNFILE_H_ 
#define RUNFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <gmpxx.h>
#include "mpfc_class.h"

class Runfile_c{
	const std::string filename;
	std::vector<std::string> lines;

	public:
		std::vector<std::vector<mpfc_class>> runs;
		std::vector<int> maxOrders;

		Runfile_c(const char* filename);
		Runfile_c(const std::string filename);
		Runfile_c(const std::vector<std::string> line);
		
		int ReadRunfile();

		int Expand();

		int ExpandBraces(const int param);

		std::tuple<mpfc_class, mpfc_class, mpfc_class> ParseBraces(std::string insideBraces);

		std::tuple<size_t, size_t> FindNthParameter(const std::string line, const int param);

		int ExpandRelativeEqns(const int param);
			
		std::tuple<mpfc_class, int> ParseRelativeEqn(std::string equation, std::string relTo);

		mpfc_class RelativeMPF(std::string firstHalf, std::string equation);

		std::string FindBaseNumber(std::string sourceString, const int paramNumber);

		int ReadMPF(mpfc_class& output, FILE* runfile);

		int ReadMaxOrder(FILE* runfile);

		int RunCompare(std::vector<mpfc_class> run1, std::vector<mpfc_class> run2);

		std::string NameOutputFile();
};

#endif
