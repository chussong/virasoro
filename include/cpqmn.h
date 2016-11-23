#ifndef CPQMN_H_
#define CPQMN_H_

#include <iostream>
#include <gmpxx.h>

template<class T>
class Cpqmn_c{
	const T* bsq, *invBsq;	
	const int numberOfMN;
	const unsigned short int maxOrder;
	const unsigned short int* mTable, *nTable;
	const int* mnLookup;

	public:
		T* hpmn;
		T* Amn;
		T* Rmn;

		Cpqmn_c(const T* bsq, const T* invBsq, const int numberOfMN, const unsigned short int maxOrder, const unsigned short int* mTable, const unsigned short int* nTable, const int* mnLookup);
		~Cpqmn_c();
		
		void FillHpmn ();

		void FillAmn();

		void FindAmn(const unsigned int m, const unsigned int n, T& temp1, T& temp2, T& prod);

		void FillRmn(const T* llsq, const T* lhsq);
};

template<class T>
Cpqmn_c<T>::Cpqmn_c(const T* bsq, const T* invBsq, const int numberOfMN, const unsigned short int maxOrder, const unsigned short int* mTable, const unsigned short int* nTable, const int* mnLookup): bsq(bsq), invBsq(invBsq), numberOfMN(numberOfMN), maxOrder(maxOrder), mTable(mTable), nTable(nTable), mnLookup(mnLookup)
{
	hpmn = new T[numberOfMN];
	Amn = new T[numberOfMN];
	Rmn = new T[numberOfMN];
}

template<class T>
Cpqmn_c<T>::~Cpqmn_c(){
	delete[] hpmn;
	delete[] Amn;
	delete[] Rmn;
//	hpmn = nullptr;
//	Amn = nullptr;
//	Rmn = nullptr;
}

template<class T>
void Cpqmn_c<T>::FillHpmn (){
	T temp;
	for(int pos = 1; pos <= numberOfMN; ++pos){
		hpmn[pos-1] = 1;
		hpmn[pos-1] -= nTable[pos-1]*nTable[pos-1];
		hpmn[pos-1] *= *bsq;
		hpmn[pos-1] /= 4;
		temp = 1;
		temp -= mTable[pos-1]*mTable[pos-1];
		temp *= *invBsq;
		temp /= 4;
		hpmn[pos-1] += temp;
		temp = 1;
		temp /= 2;
		hpmn[pos-1] += temp;
		temp = mTable[pos-1]*nTable[pos-1];
		temp /= 2;
		hpmn[pos-1] -= temp;
	}
	return;
}

template<class T>
inline void Cpqmn_c<T>::FillAmn(){
	T temp1, temp2, prod;
	for(int pos = 1; pos <= numberOfMN; ++pos){
		FindAmn(mTable[pos-1], nTable[pos-1], temp1, temp2, Amn[pos-1]);
	}
}

template<class T>
void Cpqmn_c<T>::FindAmn(const unsigned int m, const unsigned int n, T& temp1, T& temp2, T& prod){
	prod = 1;
	for(unsigned int l = 1; l <= n-1; ++l){ 		// pairs of Lml*Lm(-l)
		temp1 = m*m;
		temp1 *= *invBsq;
		temp2 = l*l;
		temp2 *= *bsq;
		temp1 -= temp2;
		prod *= temp1;
	}
	for(unsigned int k = 1; k <= m-1; ++k){ 		// pairs of Lkn*L(-k)n
		temp1 = k*k;
		temp1 *= *invBsq;
		temp2 = n*n;
		temp2 *= *bsq;
		temp2 -= temp1;
		prod *= temp2;
	}
	for(unsigned int k = 1; k < m; ++k){ 		// pairs of Lk0*L(-k)0
		temp1 = k*k;
		temp1 *= *invBsq;
		temp1 = -temp1;
		prod *= temp1;
	}
	for(unsigned int l = 1; l < n; ++l){ 		// pairs of L0l*L0(-l)
		temp1 = l*l;
		temp1 *= *bsq;
		temp1 = -temp1;
		prod *= temp1;
	}
	prod *= m*n;				// loose Lm0 and L0n
	prod = -prod;				
	for(unsigned int k = 1; k <= m-1; ++k){		
		for(unsigned int l = 1; l <= n-1; ++l){
			temp1 = k*k;
			temp1 *= *invBsq;
			temp2 = l*l;
			temp2 *= *bsq;
			temp1 += temp2;
			temp2 = temp1 + 2*k*l;					// paired L(-k)l*Lk(-l)
			temp1 -= 2*k*l;							// paired Lkl*L(-k)(-l) 
			prod *= temp1;
			prod *= temp2;
		}
	}
	if(prod == 0){
//		std::cout << "A[" << m << "," << n << "] = 1/0. This would be multiplied by P[" << m << "," << n << "] = " << this->Rmn[mnLookup[(m-1)*maxOrder + n-1]-1] << std::endl;
		return;
	}
	prod = 1/prod;
	return;
}

template<class T>
void Cpqmn_c<T>::FillRmn(const T* llsq, const T* lhsq){
	T temp1, temp2, Lsq;
	
	temp1 = *bsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	Rmn[0] = temp1*temp2;
	Rmn[0] *= *bsq;
	Rmn[0] *= *bsq;
	Rmn[0] /= 2;
	Rmn[0] = -Rmn[0];																			// R12

	temp1 = *invBsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	Rmn[1] = temp1*temp2;
	Rmn[1] *= *invBsq;
	Rmn[1] *= *invBsq;
	Rmn[1] /= 2;
	Rmn[1] = -Rmn[1];																			// R21

	temp1 = 2;
	temp1 += *invBsq;
	temp1 += *bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 += *invBsq;
	temp2 += *bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 += *invBsq;
	temp2 += *bsq;
	Rmn[3] = temp1*temp2;
	Rmn[3] *= temp2;
	temp1 = 2;
	temp1 = -temp1;
	temp1 += *invBsq;
	temp1 += *bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += *invBsq;
	temp2 += *bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += *invBsq;
	temp2 += *bsq;
	temp1 *= temp2;
	Rmn[3] *= temp1;
	Rmn[3] *= temp2;
	Rmn[3] /= 2;
	Rmn[3] = -Rmn[3];																			// R22
	
	for(int m = 1; m <= 2; ++m){
		for(int n = 3+(m%2); m*n <= maxOrder; n+=(1+m%2)){
			Rmn[mnLookup[(m-1)*maxOrder + n-1]-1] = Rmn[mnLookup[(m-1)*maxOrder + n-3]-1];
			for(int p = -m+1; p <= m-1; p+=2){
				Lsq = 2*n-1;
				Lsq -= n*n;
				Lsq *= *bsq;
				temp1 = 2*(n-1)*p;
				Lsq += temp1;
				temp1 = (*invBsq)*(p*p);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[mnLookup[(m-1)*maxOrder + n-1]-1] *= temp1;
			}
		}
	}

	for(int m = 3; m <= maxOrder-(maxOrder%2); m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){
			Rmn[mnLookup[(m-1)*maxOrder + n-1]-1] = Rmn[mnLookup[(m-3)*maxOrder + n-1]-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = 2*m-1;
				Lsq -= m*m;
				Lsq *= *invBsq;
				temp1 = 2*(m-1)*q;
				Lsq += temp1;
				temp1 = (*bsq)*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[mnLookup[(m-1)*maxOrder + n-1]-1] *= temp1;
			}
		}
		for(int n = 1; m*n+n <= maxOrder; ++n){
			Rmn[mnLookup[m*maxOrder + n-1]-1] = Rmn[mnLookup[(m-2)*maxOrder + n-1]-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = (*invBsq)*(m*m);
				Lsq = -Lsq;
				temp1 = 2*q*m;
				Lsq -= temp1;
				temp1 = (*bsq)*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[mnLookup[m*maxOrder + n-1]-1] *= temp1;
			}
		}
	}
	
	this->FillAmn();
	for(int pos = 1; pos <= numberOfMN; ++pos){
//		std::cout << "Pmn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Rmn[pos-1], 10) << std::endl;
//		std::cout << "Amn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Amn[pos-1], 10) << std::endl;
		Rmn[pos-1] *= Amn[pos-1];
//		std::cout << "Rmn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Rmn[pos-1], 10) << std::endl;
	}
	
	return;
}
#endif
