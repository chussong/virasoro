#ifndef CPQMN_H_
#define CPQMN_H_

#include <iostream>
#include <gmp.h>

class Cpqmn_t{
	const bool complex = false;
	const mpf_t* bsq, *invBsq;	
	const int numberOfMN;
	const unsigned short int maxOrder;
	const unsigned short int* mTable, *nTable;
	const int* mnLookup;

	public:
		mpf_t* hpmn;
		mpf_t* Amn;
		mpf_t* Rmn;

		Cpqmn_t(const mpf_t* bsq, const mpf_t* invBsq, const int numberOfMN, const unsigned short int maxOrder, const unsigned short int* mTable, const unsigned short int* nTable, const int* mnLookup);
		~Cpqmn_t();
		
		void FillHpmn ();

		void FillAmn();

		void FindAmn(const unsigned int m, const unsigned int n, mpf_t& temp1, mpf_t& temp2, mpf_t& prod);

		void FillRmn(const mpf_t* llsq, const mpf_t* lhsq);
};

#endif
