#ifndef VIRASORO_H_
#define VIRASORO_H_

#include <chrono>		// timers
#include <iostream>		// cout
#include <fstream>		// file output
#include <string>		// std::string
#include <vector>		// std::vector

#include "runfile.h"

int core(int argc, char** argv, const bool wolfram);
std::vector<std::string> CollectArgs(int argc, char** argv);
void ReadDefaults(const std::string filename, const bool quiet);
void CreateConfigFile(const std::string filename);
std::string ParseOptions(std::vector<std::string> &args);

std::string to_string(const mpfr::mpreal N, int digits);
std::string to_string(const std::complex<mpfr::mpreal> N, int digits, int base = 10);
#endif
