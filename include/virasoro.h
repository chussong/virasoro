#ifndef VIRASORO_H_
#define VIRASORO_H_
//#define _CRT_DISABLE_PERFCRIT_LOCKS			requires removing iostream and possibly extra static linking

#include <cstdlib>		// atoi
#include <chrono>		// timers
#include <iostream>		// cout
#include <string>
#include <fstream>		// file output
#include <stdio.h>		// no idea why I needed this
#include <gmpxx.h>		// mpf_class
#include <thread>		// std::thread

#include "cpqmn.h"
#include "hmn.h"

extern const int maxThreads;
extern const int precision;
extern mpf_class* powOverflow;

inline void ClearStructureChars(FILE* file){
	char c = fgetc(file);
	while(c == ' ' || c == ',' || c == ';') c = fgetc(file);
	ungetc(c, file);
}

int ReadRunfile(char* filename, mpf_class** &runs, int* &maxOrders);

int ReadMPF(mpf_class& output, FILE* runfile);

int ReadMaxOrder(FILE* runfile);

void SetPowOverflow(unsigned short int maxOrder);

void DebugPrintRunVector(const mpf_class* runVector, const unsigned short int maxOrder);

void FindCoefficients(const mpf_class* runVector, const unsigned short int maxOrder, const std::string runfileName);

int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ConvertInputs(mpf_class& bsq, mpf_class& invBsq, mpf_class& llsq, mpf_class& lhsq, const mpf_class& c, const mpf_class& hl, const mpf_class& hh, mpf_class& temp1, mpf_class& temp2);

void FillH(mpf_class* H, const Hmn_t* Hmn, const Cpqmn_t* Cpqmn, const mpf_class hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart);

std::string to_string(const mpf_class N, int digits);

std::string NameOutputFile(const char* runfileName);

void DisplayH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder);

void WriteH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder, const std::string runfileName);

#endif
