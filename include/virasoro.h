#include <cstdlib>		// atoi
#include <chrono>		// timers
#include <iostream>		// cout
#include <fstream>		// file output
#include <stdio.h>
#include <gmp.h>
#include <thread>

const static int maxThreads = 8;
static const int precision = 512;
static mpf_t* powOverflow;

int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void FillHpmn (mpf_t* hpmn, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN, const mpf_t bsq, const mpf_t invBsq);

void FillProdLkl(mpf_t *prodLkl, const mpf_t bsq, const mpf_t invBsq, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN);

void FindProdLkl(const mpf_t bsq, const mpf_t invBsq, const unsigned int m, const unsigned int n, mpf_t& temp1, mpf_t& temp2, mpf_t& prod);

void FillRmn(mpf_t* Rmn, const mpf_t* prodLkl, const mpf_t bsq, const mpf_t invBsq, const mpf_t llsq, const mpf_t lhsq, const int numberOfMN, const int* mnLookup, const unsigned short int maxOrder);

void FillHmn(mpf_t** Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ThreadFillHmn(mpf_t** Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const int* mnLocation, const int* mnMultiplicity, const int startingMN, const int endingMN, const int order, mpf_t& temp1, mpf_t& temp2, mpf_t& hpTemp);

void FillH(mpf_t* H, const mpf_t* const* Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const mpf_t hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void DisplayH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder, int time, const std::string unit);
