#ifndef CPQMN_H_
#define CPQMN_H_

#include <vector>
#include <iostream>
#include "access.h"

namespace virasoro {
extern mpfr::mpreal tolerance;

template<class T>
class Cpqmn_c{
	const T bsq;
	const T invBsq;
	const int maxOrder;
	Cpqmn_c() = default;

	std::vector<T> hpmn;
	std::vector<T> Amn;
	std::vector<T> Rmn;
	std::vector<std::vector<T>> Cpqmn;

	void FillAmn();
	void FindAmn(T& prod, const int m, const int n, T& temp1, T& temp2);

	public:
		Cpqmn_c(const T& bsq, const T& invBsq, const int maxOrder);
		
		static void FillHpmn(Cpqmn_c<T>& C);
		static void FillRmn(Cpqmn_c<T>& C, const T* llsq, const T* lhsq);
		static void FillCpqmn(Cpqmn_c<T>& C);

		static int CheckForDivergences(const Cpqmn_c<T>& Cpqmn, const bool showProgressBar, const bool quiet);

		const std::vector<T>& operator[](const int pos) const { return Cpqmn[pos]; }
		T CFromHp(const T& hp, const int pos) const;

		void DoneWithOrder(const int order);
};

template<class T>
Cpqmn_c<T>::Cpqmn_c(const T& bsq, const T& invBsq, const int maxOrder): 
	bsq(bsq), invBsq(invBsq), maxOrder(maxOrder)
{
	static_assert(std::is_same<T,mpfr::mpreal>::value || 
		std::is_same<T,std::complex<mpfr::mpreal>>::value, 
		"Cpqmn_c must be instantiated as Cpqmn_c<mpreal> or Cpqmn_c<complex<mpreal>>.");
	hpmn.resize(Access::TotalMN());
	Amn.resize(Access::TotalMN());
	Rmn.resize(Access::TotalMN());
	Cpqmn.resize(Access::TotalMN() - Access::MultOfMN(maxOrder));
	int rowLength = Cpqmn.size();
	int i = 0;
	for(int mn = 2; mn < maxOrder; mn += 2){
		for(int j = 0; j < Access::MultOfMN(mn); ++j){
//			std::cout << "About to resize Cpqmn[" << i << "] to " << rowLength << ".\n";
			Cpqmn[i].resize(rowLength);
			++i;
		}
		rowLength -= Access::MultOfMN(maxOrder - mn);
	}
}

template<class T>
void Cpqmn_c<T>::FillHpmn (Cpqmn_c<T>& C){
	T temp;
	for(unsigned int i = 0; i < C.hpmn.size(); ++i){
		C.hpmn[i] = 1;
		C.hpmn[i] -= Access::nAtLoc(i)*Access::nAtLoc(i);
		C.hpmn[i] *= C.bsq;
		C.hpmn[i] /= 4;
		temp = 1;
		temp -= Access::mAtLoc(i)*Access::mAtLoc(i);
		temp *= C.invBsq;
		temp /= 4;
		C.hpmn[i] += temp;
		temp = 1;
		temp /= 2;
		C.hpmn[i] += temp;
		temp = Access::mAtLoc(i)*Access::nAtLoc(i);
		temp /= 2;
		C.hpmn[i] -= temp;
	}
	return;
}

template<class T>
inline void Cpqmn_c<T>::FillAmn(){
	T temp1, temp2, prod;
	for(unsigned int i = 0; i < Amn.size(); ++i){
		FindAmn(Amn[i], Access::mAtLoc(i), Access::nAtLoc(i), temp1, temp2);
	}
}

// why the fuck am I not doing this recursively
template<class T>
void Cpqmn_c<T>::FindAmn(T& prod, const int m, const int n, T& temp1, T& temp2){
	prod = 1;
	for(int l = 1; l <= n-1; ++l){ 		// pairs of Lml*Lm(-l)
		temp1 = m*m;
		temp1 *= invBsq;
		temp2 = l*l;
		temp2 *= bsq;
		temp1 -= temp2;
		prod *= temp1;
	}
	for(int k = 1; k <= m-1; ++k){ 		// pairs of Lkn*L(-k)n
		temp1 = k*k;
		temp1 *= invBsq;
		temp2 = n*n;
		temp2 *= bsq;
		temp2 -= temp1;
		prod *= temp2;
	}
	for(int k = 1; k < m; ++k){ 		// pairs of Lk0*L(-k)0
		temp1 = k*k;
		temp1 *= invBsq;
		temp1 = -temp1;
		prod *= temp1;
	}
	for(int l = 1; l < n; ++l){ 		// pairs of L0l*L0(-l)
		temp1 = l*l;
		temp1 *= bsq;
		temp1 = -temp1;
		prod *= temp1;
	}
	prod *= m*n;				// loose Lm0 and L0n
	prod = -prod;				
	for(int k = 1; k <= m-1; ++k){		
		for(int l = 1; l <= n-1; ++l){
			temp1 = k*k;
			temp1 *= invBsq;
			temp2 = l*l;
			temp2 *= bsq;
			temp1 += temp2;
			temp2 = temp1 + 2*k*l;					// paired L(-k)l*Lk(-l)
			temp1 -= 2*k*l;							// paired Lkl*L(-k)(-l) 
			prod *= temp1;
			prod *= temp2;
		}
	}
	if(prod == 0){
		return;
	}
	prod = 1/prod;
	return;
}

template<class T>
void Cpqmn_c<T>::FillRmn(Cpqmn_c<T>& C, const T* llsq, const T* lhsq){
	T temp1, temp2, Lsq;
	
	// R_12
	temp1 = C.bsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	C.Rmn[0] = temp1*temp2;
	C.Rmn[0] *= C.bsq;
	C.Rmn[0] *= C.bsq;
	C.Rmn[0] /= 2;
	C.Rmn[0] = -C.Rmn[0];
	C.Rmn[0] <<= 8;

	// R_21
	temp1 = C.invBsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	C.Rmn[1] = temp1*temp2;
	C.Rmn[1] *= C.invBsq;
	C.Rmn[1] *= C.invBsq;
	C.Rmn[1] /= 2;
	C.Rmn[1] = -C.Rmn[1];
	C.Rmn[1] <<= 8;

	// R_22
	temp1 = 2;
	temp1 += C.invBsq;
	temp1 += C.bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 += C.invBsq;
	temp2 += C.bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 += C.invBsq;
	temp2 += C.bsq;
	C.Rmn[3] = temp1*temp2;
	C.Rmn[3] *= temp2;
	temp1 = 2;
	temp1 = -temp1;
	temp1 += C.invBsq;
	temp1 += C.bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += C.invBsq;
	temp2 += C.bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += C.invBsq;
	temp2 += C.bsq;
	temp1 *= temp2;
	C.Rmn[3] *= temp1;
	C.Rmn[3] *= temp2;
	C.Rmn[3] /= 2;
	C.Rmn[3] = -C.Rmn[3];
	C.Rmn[3] <<= 16;
	
	for(int m = 1; m <= 2; ++m){
		for(int n = 3+(m%2); m*n <= C.maxOrder; n+=(1+m%2)){
			C.Rmn[Access::PosFromMN(m,n)-1] = C.Rmn[Access::PosFromMN(m,n-2)-1];
			for(int p = -m+1; p <= m-1; p+=2){
				Lsq = 2*n-1;
				Lsq -= n*n;
				Lsq *= C.bsq;
				temp1 = 2*(n-1)*p;
				Lsq += temp1;
				temp1 = C.invBsq*(p*p);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				C.Rmn[Access::PosFromMN(m,n)-1] *= temp1;
			}
			C.Rmn[Access::PosFromMN(m,n)-1] <<= 8*m;
		}
	}

	for(int m = 3; m <= C.maxOrder-(C.maxOrder%2); m+=2){
		for(int n = 2; m*n <= C.maxOrder; n+=2){
			C.Rmn[Access::PosFromMN(m,n)-1] = C.Rmn[Access::PosFromMN(m-2,n)-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = 2*m-1;
				Lsq -= m*m;
				Lsq *= C.invBsq;
				temp1 = 2*(m-1)*q;
				Lsq += temp1;
				temp1 = C.bsq*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				C.Rmn[Access::PosFromMN(m,n)-1] *= temp1;
			}
			C.Rmn[Access::PosFromMN(m,n)-1] <<= 8*n;
		}
		for(int n = 1; m*n+n <= C.maxOrder; ++n){
			C.Rmn[Access::PosFromMN(m+1,n)-1] = C.Rmn[Access::PosFromMN(m-1,n)-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = C.invBsq*(m*m);
				Lsq = -Lsq;
				temp1 = 2*q*m;
				Lsq -= temp1;
				temp1 = C.bsq*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				C.Rmn[Access::PosFromMN(m+1,n)-1] *= temp1;
			}
			C.Rmn[Access::PosFromMN(m+1,n)-1] <<= 8*n;
		}
	}
	
	C.FillAmn();
	for(unsigned int i = 0; i < C.Rmn.size(); ++i){
		C.Rmn[i] *= C.Amn[i];
#ifdef VERBOSE_DEBUG
		std::cout << "Pmn[" << Access::mAtLoc(i) << "," << Access::nAtLoc(i) << "] = " << to_string(C.Rmn[i], 10) << std::endl;
		std::cout << "Amn[" << Access::mAtLoc(i) << "," << Access::nAtLoc(i) << "] = " << to_string(C.Amn[i], 10) << std::endl;
		std::cout << "Rmn[" << Access::mAtLoc(i) << "," << Access::nAtLoc(i) << "] = " << to_string(C.Rmn[i], 10) << std::endl;
#endif
	}
	
	return;
}

// this needs to not print stuff in __WSTP and __MATHEMATICA modes
template<class T>
int Cpqmn_c<T>::CheckForDivergences(const Cpqmn_c<T>& Cpqmn, const bool showProgressBar, const bool quiet){
	int oldMax = Cpqmn.maxOrder;
	int maxOrder = oldMax;
	bool problemIsAmn = false;
	for(unsigned int i = 0; i < Cpqmn.Rmn.size(); ++i){
		for(unsigned int j = 0; j < Cpqmn.Rmn.size(); ++j){
			if(mpfr::abs(Cpqmn.hpmn[i] + Access::mnAtLoc(i) - Cpqmn.hpmn[j]) < tolerance
					&& maxOrder+2 > std::max(Access::mnAtLoc(i), Access::mnAtLoc(j))){
				maxOrder = std::max(Access::mnAtLoc(i), Access::mnAtLoc(j))-2;
				problemIsAmn = false;
			}
		}
		if(Cpqmn.Amn[i] == 0 && maxOrder+2 > Access::mnAtLoc(i)){
			maxOrder = Access::mnAtLoc(i)-2;
			problemIsAmn = true;
		}
	}
	maxOrder = maxOrder - (maxOrder%2);
	if(maxOrder < oldMax){
		if(showProgressBar) std::cout << "\r";
		if(!quiet){
			if(maxOrder > 2) std::cout << "Stopping this run at order " 
				<< maxOrder << " because ";
			if(maxOrder <= 2) std::cout << "Skipping this run because ";

			if(!problemIsAmn) std::cout << "the coefficients diverge ";
			if(problemIsAmn) std::cout << "Amn diverges ";

			if(maxOrder > 2) std::cout << "above this." << std::endl;
			if(maxOrder <= 2) std::cout << "immediately." << std::endl;
		}
	}
	return maxOrder;
}

template<class T>
void Cpqmn_c<T>::FillCpqmn(Cpqmn_c<T>& C){
	T temp;
	for(unsigned int i = 0; i < C.Cpqmn.size(); ++i){
		for(unsigned int j = 0; j < C.Cpqmn[i].size(); ++j){
//			std::cout << "Filling Cpqmn[" << i << "," << j << "]." << std::endl;
			temp = C.hpmn[i] + Access::mAtLoc(i)*Access::nAtLoc(i);
			C.Cpqmn[i][j] = C.Rmn[j];
			temp -= C.hpmn[j];
			C.Cpqmn[i][j] /= temp;
		}
	}
	return;
}

template<class T>
T Cpqmn_c<T>::CFromHp(const T& hp, const int pos) const {
	T ret = hp - hpmn[pos-1];
	ret = Rmn[pos-1]/ret; // my Rmn includes the factor 16^mn
	return ret;
}

template<class T>
void Cpqmn_c<T>::DoneWithOrder(const int order){
	Cpqmn.resize(Cpqmn.size() - Access::MultOfMN(maxOrder - order));
	Cpqmn.shrink_to_fit();
}
} // namespace virasoro
#endif
