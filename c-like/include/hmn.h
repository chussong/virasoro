#ifndef HMN_H_
#define HMN_H_

#include <thread>
#include <gmp.h>
#include "cpqmn.h"

extern const int maxThreads;
extern const int precision;
extern mpf_t* powOverflow;

class Hmn_t{
	const int numberOfMN;
	const unsigned short int maxOrder;
	const int* mnLocation, *mnMultiplicity;
	const int* mnLookup;
	int* mnAtThisOrder;
	
	public:
		mpf_t** Hmn;
		Cpqmn_t* Cpqmn;
		
		Hmn_t(Cpqmn_t* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup);
		~Hmn_t();
		
		inline static void StartThread(Hmn_t* Hmn, const int startingMN, const int endingMN, const int order, mpf_t& temp1, mpf_t& temp2, mpf_t& hpTemp){
			Hmn->ThreadFillHmn(startingMN, endingMN, order, temp1, temp2, hpTemp);
		}
		
		void FillHmn();

		void ThreadFillHmn(const int startingMN, const int endingMN, const int order, mpf_t& temp1, mpf_t& temp2, mpf_t& hpTemp);
};
#endif
