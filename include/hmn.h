#ifndef HMN_H_
#define HMN_H_

#include <thread>
#include <gmpxx.h>
#include "cpqmn.h"

extern int maxThreads;
extern int precision;
extern bool showProgressBar;

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
		
		inline static void StartThread(Hmn_c<T>* Hmn, const int startingMN, const int endingMN, const int order, T& temp){
			Hmn->ThreadFillHmn(startingMN, endingMN, order, temp);
		}
		
		void FillHmn();

		void ThreadFillHmn(const int startingMN, const int endingMN, const int order, T& temp);
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

inline void DrawProgressBar(const float progress){
	int barWidth = 70;
	std::cout << "[";
	int progPos = barWidth*progress;
	for(int i = 1; i <= barWidth; ++i){
		if(i <= progPos) std::cout << "-";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();
}

template<class T>
void Hmn_c<T>::FillHmn(){	
	std::thread* thread = new std::thread[maxThreads];
	T* temp = new T[maxThreads];
	int numThreads;
	float progress = 0.0;
	float totalComputations = maxOrder*maxOrder*maxOrder/6.0f; 
	if(showProgressBar){
		for(int mn = 2; mn < maxOrder; mn+=2){
			totalComputations += 0.5f*(maxOrder-mn)*mnMultiplicity[mn-1];
		}
		std::cout << "\r";
		DrawProgressBar(progress);
	}
	for(int order = 2; order < maxOrder; order+=2){	
		numThreads = std::min(maxThreads,(maxOrder-order)/2);
		thread[0] = std::thread(Hmn_c::StartThread, this, 2, (maxOrder-order)/numThreads + (maxOrder-order)%numThreads, order, std::ref(temp[0]));
		for(int i=2; i<=numThreads; ++i){	
			thread[i-1] = std::thread(Hmn_c::StartThread, this, (i-1)*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads + 1, i*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads, order, std::ref(temp[i-1]));
		}
		for(int i=1; i<= numThreads; ++i){
			thread[i-1].join();
		}
		if(showProgressBar){
//			progress += (float)2/maxOrder;
			progress += order*(maxOrder-order)/totalComputations;
			DrawProgressBar(progress);
		}
	}
	if(showProgressBar){
		progress = 1.0f;
		DrawProgressBar(progress);
		std::cout << std::endl;
	}
	delete[] temp;
	delete[] thread;
	return;
}

// change this to only use pos
template<class T>
void Hmn_c<T>::ThreadFillHmn(const int startingMN, const int endingMN, const int order, T& temp){
	for(int mn = startingMN + startingMN%2; mn <= endingMN; mn+=2){
		for(int pos = mnLocation[mn-1]; pos <= mnLocation[mn-1] + mnMultiplicity[mn-1] - 1; ++pos){
			for(int power = 2; power <= order; power+=2){
				for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
					temp = Cpqmn->Cpqmn[pos-1][scanPos-1];
					temp *= Hmn[(order-power)/2][scanPos-1];
					Hmn[order/2][pos-1] += temp;
				}
			}
		}
	}
	return;
}
#endif
