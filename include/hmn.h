#ifndef HMN_H_
#define HMN_H_

#include <thread>
#include <gmpxx.h>
#include "cpqmn.h"

extern int maxThreads;
extern int precision;
extern mpf_class* powOverflow;

template <class T> // to be used with either mpf_class or mpfc_class
class Hmn_c{
	const int numberOfMN;
	const unsigned short int maxOrder;
	const int* mnLocation, *mnMultiplicity;
	const int* mnLookup;
	int* mnAtThisOrder;
	
	public:
		T** Hmn;
		Cpqmn_c<T>* Cpqmn;
		
		Hmn_c(Cpqmn_c<T>* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup);
		~Hmn_c();
		
		inline static void StartThread(Hmn_c<T>* Hmn, const int startingMN, const int endingMN, const int order, T& temp1, T& temp2, T& hpTemp){
			Hmn->ThreadFillHmn(startingMN, endingMN, order, temp1, temp2, hpTemp);
		}
		
		void FillHmn();

		void ThreadFillHmn(const int startingMN, const int endingMN, const int order, T& temp1, T& temp2, T& hpTemp);
};

template<class T>
Hmn_c<T>::Hmn_c(Cpqmn_c<T>* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup): numberOfMN(numberOfMN), maxOrder(maxOrder), mnLocation(mnLocation), mnMultiplicity(mnMultiplicity), mnLookup(mnLookup), Cpqmn(Cpqmn)
{
	Hmn = new T*[maxOrder/2];
	Hmn[0] = new T[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos) Hmn[0][pos-1] = 1;
	mnAtThisOrder = new int[maxOrder/2];
	for(int order = 2; order <= maxOrder-2; order+=2){
		mnAtThisOrder[order/2-1] = 0;
		for(int mn = 2; mn <= maxOrder-order; mn+=2){
			mnAtThisOrder[order/2-1] += mnMultiplicity[mn-1];
		}
		Hmn[order/2] = new T[mnAtThisOrder[order/2-1]];
	}
}

template<class T>
Hmn_c<T>::~Hmn_c(){
	delete[] Hmn[0];
	for(int order = 2; order <= maxOrder-2; order +=2){
		delete[] Hmn[order/2];
	}
	delete[] Hmn;
	delete[] mnAtThisOrder;
}

template<class T>
void Hmn_c<T>::FillHmn(){	
	std::thread* thread = new std::thread[maxThreads];
//	T temp1[maxThreads], temp2[maxThreads], hpTemp[maxThreads];
	T* temp1 = new T[maxThreads];
	T* temp2 = new T[maxThreads];
	T* hpTemp = new T[maxThreads];
	int numThreads;
	for(int order = 2; order < maxOrder; order+=2){	
		numThreads = std::min(maxThreads,(maxOrder-order)/2);
		thread[0] = std::thread(Hmn_c::StartThread, this, 2, (maxOrder-order)/numThreads + (maxOrder-order)%numThreads, order, std::ref(temp1[0]), std::ref(temp2[0]), std::ref(hpTemp[0]));
		for(int i=2; i<=numThreads; ++i){	
			thread[i-1] = std::thread(Hmn_c::StartThread, this, (i-1)*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads + 1, i*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads, order, std::ref(temp1[i-1]), std::ref(temp2[i-1]), std::ref(hpTemp[i-1]));
		}
		for(int i=1; i<= numThreads; ++i){
			thread[i-1].join();
		}
	}
	delete[] temp1;
	delete[] temp2;
	delete[] hpTemp;
	delete[] thread;
	return;
}

template<class T>
void Hmn_c<T>::ThreadFillHmn(const int startingMN, const int endingMN, const int order, T& temp1, T& temp2, T& hpTemp){
	for(int mn = startingMN + startingMN%2; mn <= endingMN; mn+=2){
		for(int pos = mnLocation[mn-1]; pos <= mnLocation[mn-1] + mnMultiplicity[mn-1] - 1; ++pos){
			hpTemp = Cpqmn->hpmn[pos-1] + mn;
			for(int power = 2; power <= order; power+=2){
				for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
					temp1 = Cpqmn->Rmn[scanPos-1];
					temp1 *= Hmn[(order-power)/2][scanPos-1];
/*					temp1 *= powOverflow[power/256];
					temp2 = 16;
					mpf_pow_ui(temp2.get_mpf_t(), temp2.get_mpf_t(), power%256);
					temp1 *= temp2;*/
					temp1 <<= 4*power;		// fast multiplication by 2^(4*power)
					temp2 = hpTemp - Cpqmn->hpmn[scanPos-1];
/*					if(temp2 == 0){
						std::cout << "At order " << order << ", mn = " << mn <<", dividing numerator " << temp1 << " by " << temp2 << "." << std::endl;
						continue;
					}*/
					temp1 /= temp2;
					Hmn[order/2][pos-1] += temp1;
				}
			}
		}
	}
	return;
}
#endif
