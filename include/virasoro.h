#ifndef VIRASORO_H_
#define VIRASORO_H_
//#define _CRT_DISABLE_PERFCRIT_LOCKS			requires removing iostream and possibly extra static linking

#include <cstdlib>		// atoi
#include <chrono>		// timers
#include <iostream>		// cout
#include <string>		// std::string
#include <vector>		// std::vector
#include <fstream>		// file output
#include <stdio.h>		// fgetc
#include <gmpxx.h>		// mpf_class
#include <thread>		// std::thread
#include <tuple>		// std::tuple

#include "runfile.h"

std::vector<std::string> CollectArgs(int argc, char** argv);

void ReadDefaults(std::string filename);

void CreateConfigFile(std::string filename);

std::string ParseOptions(std::vector<std::string> &args);

#endif
