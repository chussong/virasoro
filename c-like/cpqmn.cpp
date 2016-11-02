#include "cpqmn.h"

Cpqmn_t::Cpqmn_t(const mpf_t* bsq, const mpf_t* invBsq, const int numberOfMN, const unsigned short int maxOrder, const unsigned short int* mTable, const unsigned short int* nTable, const int* mnLookup): bsq(bsq), invBsq(invBsq), numberOfMN(numberOfMN), maxOrder(maxOrder), mTable(mTable), nTable(nTable), mnLookup(mnLookup)
{
	hpmn = new mpf_t[numberOfMN];
	Amn = new mpf_t[numberOfMN];
	Rmn = new mpf_t[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_inits(hpmn[pos-1], Amn[pos-1], Rmn[pos-1], NULL);
	}
}

Cpqmn_t::~Cpqmn_t(){
	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_clears(hpmn[pos-1], Amn[pos-1], Rmn[pos-1], NULL);
	}
	delete[] hpmn;
	delete[] Amn;
	delete[] Rmn;
	hpmn = nullptr;
	Amn = nullptr;
	Rmn = nullptr;
}

void Cpqmn_t::FillHpmn (){
	mpf_t temp;
	mpf_init(temp);
	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_init(hpmn[pos-1]);
		mpf_set_ui(hpmn[pos-1], 1);
		mpf_sub_ui(hpmn[pos-1], hpmn[pos-1], nTable[pos-1]*nTable[pos-1]);
		mpf_mul(hpmn[pos-1], hpmn[pos-1], *bsq);
		mpf_div_ui(hpmn[pos-1], hpmn[pos-1], 4);
		mpf_set_ui(temp, 1);
		mpf_sub_ui(temp, temp, mTable[pos-1]*mTable[pos-1]);
		mpf_mul(temp, temp, *invBsq);
		mpf_div_ui(temp, temp, 4);
		mpf_add(hpmn[pos-1], hpmn[pos-1], temp);
		mpf_set_ui(temp, 1);
		mpf_div_ui(temp, temp, 2);
		mpf_add(hpmn[pos-1], hpmn[pos-1], temp);
		mpf_set_ui(temp, mTable[pos-1]*nTable[pos-1]);
		mpf_div_ui(temp, temp, 2);
		mpf_sub(hpmn[pos-1], hpmn[pos-1], temp);
	}
	mpf_clear(temp);
	return;
}

void Cpqmn_t::FillAmn(){
	mpf_t temp1, temp2, prod;
	mpf_inits(temp1, temp2, prod, NULL);
	for(int pos = 1; pos <= numberOfMN; ++pos){
		FindAmn(mTable[pos-1], nTable[pos-1], temp1, temp2, prod);
		mpf_set(Amn[pos-1], prod);
	}
	mpf_clears(temp1, temp2, prod, NULL);
}

void Cpqmn_t::FindAmn(const unsigned int m, const unsigned int n, mpf_t& temp1, mpf_t& temp2, mpf_t& prod){
	mpf_set_ui(prod,1);
	for(unsigned int l = 1; l <= n-1; ++l){ 		// pairs of Lml*Lm(-l)
		mpf_set_ui(temp1, m*m);
		mpf_mul(temp1, temp1, *invBsq);
		mpf_set_ui(temp2, l*l);
		mpf_mul(temp2, temp2, *bsq);
		mpf_sub(temp1, temp1, temp2);
		mpf_mul(prod, prod, temp1);
	}
	for(unsigned int k = 1; k <= m-1; ++k){ 		// pairs of Lkn*L(-k)n
		mpf_set_ui(temp1, k*k);
		mpf_mul(temp1, temp1, *invBsq);
		mpf_set_ui(temp2, n*n);
		mpf_mul(temp2, temp2, *bsq);
		mpf_sub(temp1, temp2, temp1);
		mpf_mul(prod, prod, temp1);
	}
	for(unsigned int k = 1; k < m; ++k){ 		// pairs of Lk0*L(-k)0
		mpf_set_si(temp1, k*k);
		mpf_mul(temp1, temp1, *invBsq);
		mpf_neg(temp1, temp1);
		mpf_mul(prod, prod, temp1);		
	}
	for(unsigned int l = 1; l < n; ++l){ 		// pairs of L0l*L0(-l)
		mpf_set_ui(temp1, l*l);
		mpf_mul(temp1, temp1, *bsq);
		mpf_neg(temp1, temp1);
		mpf_mul(prod, prod, temp1);	
	}
	mpf_mul_ui(prod, prod, m*n);
	mpf_neg(prod, prod);				// loose Lm0 and L0n
	for(unsigned int k = 1; k <= m-1; ++k){		
		for(unsigned int l = 1; l <= n-1; ++l){
			mpf_set_ui(temp1, k*k);
			mpf_mul(temp1, temp1, *invBsq);
			mpf_set_ui(temp2, l*l);
			mpf_mul(temp2, temp2, *bsq);
			mpf_add(temp1, temp1, temp2);
			mpf_add_ui(temp2, temp1, 2*k*l);		// paired L(-k)l*Lk(-l)
			mpf_sub_ui(temp1, temp1, 2*k*l);		// paired Lkl*L(-k)(-l) 
			mpf_mul(prod, prod, temp1);
			mpf_mul(prod, prod, temp2);					
		}
	}
	mpf_ui_div(prod, 1, prod);
	return;
}

void Cpqmn_t::FillRmn(const mpf_t* llsq, const mpf_t* lhsq){
	mpf_t temp1, temp2, Lsq;
	mpf_inits(temp1, temp2, Lsq, NULL);
	
	mpf_div_ui(temp1, *bsq, 16);
	mpf_set(temp2, temp1);
	mpf_add(temp1, temp1, *llsq);
	mpf_add(temp2, temp2, *lhsq);
	mpf_mul(Rmn[0], temp1, temp2);
	mpf_mul(Rmn[0], Rmn[0], *bsq);
	mpf_mul(Rmn[0], Rmn[0], *bsq);
	mpf_div_ui(Rmn[0], Rmn[0], 2);
	mpf_neg(Rmn[0], Rmn[0]);																	// R12

	mpf_div_ui(temp1, *invBsq, 16);
	mpf_set(temp2, temp1);
	mpf_add(temp1, temp1, *llsq);
	mpf_add(temp2, temp2, *lhsq);
	mpf_mul(Rmn[1], temp1, temp2);
	mpf_mul(Rmn[1], Rmn[1], *invBsq);
	mpf_mul(Rmn[1], Rmn[1], *invBsq);
	mpf_div_ui(Rmn[1], Rmn[1], 2);
	mpf_neg(Rmn[1], Rmn[1]);																	// R21

	mpf_set_ui(temp1, 2);
	mpf_add(temp1, temp1, *invBsq);
	mpf_add(temp1, temp1, *bsq);
	mpf_div_ui(temp1, temp1, 16);
	mpf_add(temp1, temp1, *llsq);
	mpf_set_ui(temp2, 2);
	mpf_add(temp2, temp2, *invBsq);
	mpf_add(temp2, temp2, *bsq);
	mpf_div_ui(temp2, temp2, 16);
	mpf_add(temp2, temp2, *lhsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_set_ui(temp2, 2);
	mpf_add(temp2, temp2, *invBsq);
	mpf_add(temp2, temp2, *bsq);
	mpf_mul(Rmn[3], temp1, temp2);
	mpf_mul(Rmn[3], Rmn[3], temp2);
	mpf_set_ui(temp1, 2);
	mpf_neg(temp1, temp1);
	mpf_add(temp1, temp1, *invBsq);
	mpf_add(temp1, temp1, *bsq);
	mpf_div_ui(temp1, temp1, 16);
	mpf_add(temp1, temp1, *llsq);
	mpf_set_ui(temp2, 2);
	mpf_neg(temp2, temp2);
	mpf_add(temp2, temp2, *invBsq);
	mpf_add(temp2, temp2, *bsq);
	mpf_div_ui(temp2, temp2, 16);
	mpf_add(temp2, temp2, *lhsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_set_ui(temp2, 2);
	mpf_neg(temp2, temp2);
	mpf_add(temp2, temp2, *invBsq);
	mpf_add(temp2, temp2, *bsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_mul(Rmn[3], Rmn[3], temp1);
	mpf_mul(Rmn[3], Rmn[3], temp2);
	mpf_div_ui(Rmn[3], Rmn[3], 2);
	mpf_neg(Rmn[3], Rmn[3]);											// R22
	
	for(int m = 1; m <= 2; ++m){
		for(int n = 3+(m%2); m*n <= maxOrder; n+=(1+m%2)){
			mpf_set(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-3]-1]);
			for(int p = -m+1; p <= m-1; p+=2){
				mpf_set_ui(Lsq, 2*n-1);
				mpf_sub_ui(Lsq, Lsq, n*n);
				mpf_mul(Lsq, Lsq, *bsq);
				mpf_set_si(temp1, 2*(n-1)*p);
				mpf_add(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, *invBsq, p*p);
				mpf_sub(Lsq, Lsq, temp1);
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, *llsq);
				mpf_sub(temp2, temp2, *lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], temp1);
			}
		}
	}

	for(int m = 3; m <= maxOrder-(maxOrder%2); m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){
			mpf_set(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-3)*maxOrder + n-1]-1]);
			for(int q = -n+1; q <= n-1; q+=2){
				mpf_set_ui(Lsq, 2*m-1);
				mpf_sub_ui(Lsq, Lsq, m*m);
				mpf_mul(Lsq, Lsq, *invBsq);
				mpf_set_si(temp1, 2*(m-1)*q);
				mpf_add(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, *bsq, q*q);
				mpf_sub(Lsq, Lsq, temp1);			
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, *llsq);
				mpf_sub(temp2, temp2, *lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], temp1);
			}
		}
		for(int n = 1; m*n+n <= maxOrder; ++n){
			mpf_set(Rmn[mnLookup[m*maxOrder + n-1]-1], Rmn[mnLookup[(m-2)*maxOrder + n-1]-1]);
			for(int q = -n+1; q <= n-1; q+=2){
				mpf_mul_ui(Lsq, *invBsq, m*m);
				mpf_neg(Lsq, Lsq);
				mpf_set_si(temp1, 2*q*m);
				mpf_sub(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, *bsq, q*q);
				mpf_sub(Lsq, Lsq, temp1);
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, *llsq);
				mpf_sub(temp2, temp2, *lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[m*maxOrder + n-1]-1], Rmn[mnLookup[m*maxOrder + n-1]-1], temp1);
			}
		}
	}
	
	FILE* devnull = fopen("/dev/null", "w");
	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_mul(Rmn[pos-1], Rmn[pos-1], Amn[pos-1]);
		mpf_out_str(devnull, 10, 10, Rmn[pos-1]);		// it segfaults without this. I have no idea why.
	}
	fclose(devnull);
	
	mpf_clears(temp1, temp2, Lsq);
	return;
}
