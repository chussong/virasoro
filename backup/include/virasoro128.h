#include <cstdlib>		// atoi
#include <cmath>		// sqrt
#include <chrono>		// timers
#include <iostream>		// cout
#include <fstream>		// file output
#include <quadmath.h>	// __float128
#include <stdio.h>
#include <gmp.h>
#include <thread>

const static int maxThreads = 8;
static __float128* powOverflow;

int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);

__float128 CToB ( __float128 c);

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void FillHpmn (__float128 *hpmn, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN, const __float128 bsq, const __float128 invBsq);

void FillProdLkl(__float128 *prodLkl, const __float128 bsq, const __float128 invBsq, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN);

__float128 FindProdLkl(__float128 bsq, __float128 invBsq, int m, int n);

void FillRmn(__float128 *Rmn, const __float128 *prodLkl, const __float128 bsq, const __float128 invBsq, const __float128 llsq, const __float128 lhsq, const int numberOfMN, const int* mnLookup, const unsigned short int maxOrder);

void FillHmn(__float128** Hmn, const __float128* const* Cpqmn, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ThreadFillHmn(__float128** Hmn, const __float128* const* Cpqmn, const int* mnLocation, const int* mnMultiplicity, const int startingMN, const int endingMN, const int order);

void FillH(__float128* H, const __float128* const* Hmn, const __float128* Rmn, const __float128* hpmn, const __float128 hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void DisplayH(const __float128* H, const __float128 c, const __float128 hl, const __float128 hh, const __float128 hp, const unsigned short int maxOrder, int time, const std::string unit);

inline __float128 HmnTerm(const __float128* const* Hmn, const __float128* const* Cpqmn, const int pqPos, const int* mnLocation, const int* mnMultiplicity, const int order){
	__float128 term = 0;
	for(int power = 2; power <= order; power+=2){
		for(int mnPos = mnLocation[power/2-1]; mnPos <= mnLocation[power/2-1] + mnMultiplicity[power/2-1] - 1; ++mnPos){
			term += powOverflow[power/256]*(__float128)pow(16.0,power%256)*Cpqmn[pqPos-1][mnPos-1]*Hmn[(order-power)/2][mnPos-1];
		}
	}
	return term;
}
