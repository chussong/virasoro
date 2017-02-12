#ifndef CPQMN_H_
#define CPQMN_H_

#include <vector>
#include <iostream>
#include "access.h"

// namespace virasoro {
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

	public:
		Cpqmn_c(const T& bsq, const T& invBsq, const int maxOrder);
		
		void FillHpmn ();
		void FillAmn();
		void FindAmn(T& prod, const int m, const int n, T& temp1, T& temp2);
		void FillRmn(const T* llsq, const T* lhsq);

		static int CheckForDivergences(const Cpqmn_c<T>& Cpqmn, const bool showProgressBar, const bool quiet);

		void FillCpqmn();

		const std::vector<T>& operator[](const int pos) const { return Cpqmn[pos]; }
		T CFromHp(const T& hp, const int pos) const;
};

template<class T>
Cpqmn_c<T>::Cpqmn_c(const T& bsq, const T& invBsq, const int maxOrder): bsq(bsq), invBsq(invBsq), maxOrder(maxOrder)
{
	static_assert(std::is_same<T,mpfr::mpreal>::value || std::is_same<T,std::complex<mpfr::mpreal>>::value, "Cpqmn_c must be instantiated as Cpqmn_c<mpreal> or Cpqmn_c<complex<mpreal>>.");
	hpmn.resize(Access::TotalMN());
	Amn.resize(Access::TotalMN());
	Rmn.resize(Access::TotalMN());
	Cpqmn.resize(Access::TotalMN());
	for(int i = 0; i < Access::TotalMN(); ++i) Cpqmn[i].resize(Access::TotalMN());
}

template<class T>
void Cpqmn_c<T>::FillHpmn (){
	T temp;
	for(unsigned int i = 0; i < hpmn.size(); ++i){
		hpmn[i] = 1;
		hpmn[i] -= Access::nAtLoc(i)*Access::nAtLoc(i);
		hpmn[i] *= bsq;
		hpmn[i] /= 4;
		temp = 1;
		temp -= Access::mAtLoc(i)*Access::mAtLoc(i);
		temp *= invBsq;
		temp /= 4;
		hpmn[i] += temp;
		temp = 1;
		temp /= 2;
		hpmn[i] += temp;
		temp = Access::mAtLoc(i)*Access::nAtLoc(i);
		temp /= 2;
		hpmn[i] -= temp;
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
void Cpqmn_c<T>::FillRmn(const T* llsq, const T* lhsq){
	T temp1, temp2, Lsq;
	
	// R_12
	temp1 = bsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	Rmn[0] = temp1*temp2;
	Rmn[0] *= bsq;
	Rmn[0] *= bsq;
	Rmn[0] /= 2;
	Rmn[0] = -Rmn[0];
	Rmn[0] <<= 8;

	// R_21
	temp1 = invBsq/16;
	temp2 = temp1;
	temp1 += *llsq;
	temp2 += *lhsq;
	Rmn[1] = temp1*temp2;
	Rmn[1] *= invBsq;
	Rmn[1] *= invBsq;
	Rmn[1] /= 2;
	Rmn[1] = -Rmn[1];
	Rmn[1] <<= 8;

	// R_22
	temp1 = 2;
	temp1 += invBsq;
	temp1 += bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 += invBsq;
	temp2 += bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 += invBsq;
	temp2 += bsq;
	Rmn[3] = temp1*temp2;
	Rmn[3] *= temp2;
	temp1 = 2;
	temp1 = -temp1;
	temp1 += invBsq;
	temp1 += bsq;
	temp1 /= 16;
	temp1 += *llsq;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += invBsq;
	temp2 += bsq;
	temp2 /= 16;
	temp2 += *lhsq;
	temp1 *= temp2;
	temp2 = 2;
	temp2 = -temp2;
	temp2 += invBsq;
	temp2 += bsq;
	temp1 *= temp2;
	Rmn[3] *= temp1;
	Rmn[3] *= temp2;
	Rmn[3] /= 2;
	Rmn[3] = -Rmn[3];
	Rmn[3] <<= 16;
	
	for(int m = 1; m <= 2; ++m){
		for(int n = 3+(m%2); m*n <= maxOrder; n+=(1+m%2)){
			Rmn[Access::PosFromMN(m,n)-1] = Rmn[Access::PosFromMN(m,n-2)-1];
			for(int p = -m+1; p <= m-1; p+=2){
				Lsq = 2*n-1;
				Lsq -= n*n;
				Lsq *= bsq;
				temp1 = 2*(n-1)*p;
				Lsq += temp1;
				temp1 = invBsq*(p*p);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[Access::PosFromMN(m,n)-1] *= temp1;
			}
			Rmn[Access::PosFromMN(m,n)-1] <<= 8*m;
		}
	}

	for(int m = 3; m <= maxOrder-(maxOrder%2); m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){
			Rmn[Access::PosFromMN(m,n)-1] = Rmn[Access::PosFromMN(m-2,n)-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = 2*m-1;
				Lsq -= m*m;
				Lsq *= invBsq;
				temp1 = 2*(m-1)*q;
				Lsq += temp1;
				temp1 = bsq*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[Access::PosFromMN(m,n)-1] *= temp1;
			}
			Rmn[Access::PosFromMN(m,n)-1] <<= 8*n;
		}
		for(int n = 1; m*n+n <= maxOrder; ++n){
			Rmn[Access::PosFromMN(m+1,n)-1] = Rmn[Access::PosFromMN(m-1,n)-1];
			for(int q = -n+1; q <= n-1; q+=2){
				Lsq = invBsq*(m*m);
				Lsq = -Lsq;
				temp1 = 2*q*m;
				Lsq -= temp1;
				temp1 = bsq*(q*q);
				Lsq -= temp1;
				temp1 = Lsq/16;
				temp2 = temp1;
				temp1 -= *llsq;
				temp2 -= *lhsq;
				temp1 *= temp2;
				temp1 *= Lsq;
				temp1 *= Lsq;
				Rmn[Access::PosFromMN(m+1,n)-1] *= temp1;
			}
			Rmn[Access::PosFromMN(m+1,n)-1] <<= 8*n;
		}
	}
	
	this->FillAmn();
	for(unsigned int i = 0; i < Rmn.size(); ++i){
		Rmn[i] *= Amn[i];
#ifdef VERBOSE_DEBUG
		std::cout << "Pmn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Rmn[pos-1], 10) << std::endl;
		std::cout << "Amn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Amn[pos-1], 10) << std::endl;
		std::cout << "Rmn[" << mTable[pos-1] << "," << nTable[pos-1] << "] = " << to_string(Rmn[pos-1], 10) << std::endl;
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
			if(maxOrder > 2) std::cout << "Stopping this run at order " << maxOrder << " because ";
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
void Cpqmn_c<T>::FillCpqmn(){
	T temp;
	for(unsigned int i = 0; i < Cpqmn.size(); ++i){
		for(unsigned int j = 0; j < Cpqmn.size(); ++j){
			temp = hpmn[i] + Access::mAtLoc(i)*Access::nAtLoc(i);
			Cpqmn[i][j] = Rmn[j];
			temp -= hpmn[j];
			Cpqmn[i][j] /= temp;
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
//} // namespace virasoro
#endif
