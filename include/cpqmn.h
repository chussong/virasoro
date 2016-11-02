#ifndef CPQMN_H_
#define CPQMN_H_

#include <iostream>
#include <gmpxx.h>

class Cpqmn_t{
	const bool complex = false;
	const mpf_class* bsq, *invBsq;	
	const int numberOfMN;
	const unsigned short int maxOrder;
	const unsigned short int* mTable, *nTable;
	const int* mnLookup;

	public:
		mpf_class* hpmn;
		mpf_class* Amn;
		mpf_class* Rmn;

		Cpqmn_t(const mpf_class* bsq, const mpf_class* invBsq, const int numberOfMN, const unsigned short int maxOrder, const unsigned short int* mTable, const unsigned short int* nTable, const int* mnLookup);
		~Cpqmn_t();
		
		void FillHpmn ();

		void FillAmn();

		void FindAmn(const unsigned int m, const unsigned int n, mpf_class& temp1, mpf_class& temp2, mpf_class& prod);

		void FillRmn(const mpf_class* llsq, const mpf_class* lhsq);
};

#endif
