#ifndef HMN_H_
#define HMN_H_

#include <thread>
#include <gmpxx.h>
#include "cpqmn.h"

extern int maxThreads;
extern int precision;
extern mpf_class* powOverflow;

class Hmn_t{
	const int numberOfMN;
	const unsigned short int maxOrder;
	const int* mnLocation, *mnMultiplicity;
	const int* mnLookup;
	int* mnAtThisOrder;
	
	public:
		mpf_class** Hmn;
		Cpqmn_t* Cpqmn;
		
		Hmn_t(Cpqmn_t* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup);
		~Hmn_t();
		
		inline static void StartThread(Hmn_t* Hmn, const int startingMN, const int endingMN, const int order, mpf_class& temp1, mpf_class& temp2, mpf_class& hpTemp){
			Hmn->ThreadFillHmn(startingMN, endingMN, order, temp1, temp2, hpTemp);
		}
		
		void FillHmn();

		void ThreadFillHmn(const int startingMN, const int endingMN, const int order, mpf_class& temp1, mpf_class& temp2, mpf_class& hpTemp);
};
#endif
