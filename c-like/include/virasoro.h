#ifndef VIRASORO_H_
#define VIRASORO_H_
//#define _CRT_DISABLE_PERFCRIT_LOCKS			requires removing iostream and possibly extra static linking

#include <cstdlib>		// atoi
#include <chrono>		// timers
#include <iostream>		// cout
#include <fstream>		// file output
#include <stdio.h>		// no idea why I needed this
#include <gmp.h>		// mpf_t
#include <thread>		// std::thread

#include "cpqmn.h"
#include "hmn.h"

extern const int maxThreads;
extern const int precision;
extern mpf_t* powOverflow;

int ReadRunfile(char* filename, mpf_t** &runs, int* &maxOrders);

int ReadMPF(mpf_t& output, FILE* runfile);

int ReadMaxOrder(FILE* runfile);

void SetPowOverflow(unsigned short int maxOrder);

void DebugPrintRunVector(const mpf_t* runVector, const unsigned short int maxOrder);

void FindCoefficients(const mpf_t* runVector, const unsigned short int maxOrder, const bool multirun);

int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ConvertInputs(mpf_t& bsq, mpf_t& invBsq, mpf_t& llsq, mpf_t& lhsq, const mpf_t& c, const mpf_t& hl, const mpf_t& hh, mpf_t& temp1, mpf_t& temp2);

void FillH(mpf_t* H, const Hmn_t* Hmn, const Cpqmn_t* Cpqmn, const mpf_t hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart);

std::string to_string(const mpf_t N, int digits);

void DisplayH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder);

void WriteH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder, const bool multirun);

#endif
