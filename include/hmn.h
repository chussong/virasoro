#ifndef HMN_H_
#define HMN_H_

#include <thread>
//#include <unistd.h>
#include "cpqmn.h"
#include "access.h"

namespace virasoro {
extern int maxThreads;
extern int precision;
extern bool showProgressBar;

// the compiler should force you to instantiate this only with T=mpreal or T=complex<mpreal>
template <class T> 
class Hmn_c{
	std::vector<std::vector<T>> Diag; // diagonals consume less RAM over time than Hmn proper
	Cpqmn_c<T>* Cpqmn;
	const std::vector<T> hp;
	std::vector<int> rowLengths;
//	bool swapping;
	
	Hmn_c() = default; // default constructor is private so no one can use it

	void ThreadFillHmn(const int startingPos, const int endingPos, const int order, T& temp1);
	static void FillH(Hmn_c<T>& Hmn, const int order, T temp);

	public:
		const bool HIsComplex;
		std::vector<std::vector<mpfr::mpreal>> realH;
		std::vector<std::vector<std::complex<mpfr::mpreal>>> complexH;
		Hmn_c(Cpqmn_c<T>* Cpqmn, std::vector<T> hp, const int maxOrder, const bool complexArgs);
		
/*		inline static long GetAvailableMemory(){
			long pages = sysconf(_SC_AVPHYS_PAGES);
			long page_size = sysconf(_SC_PAGESIZE);
			return pages*page_size;
		}*/

		std::vector<T> operator[] (int orderBy2) { return Diag[orderBy2]; } // is this needed?
		inline T& at(int order, int loc) { return Diag[order/2][loc]; } // not sure we need this either
		inline T& at(int order, int m, int n) { return Diag[order/2 + m*n/2 - 1][Access::PosFromMN(m,n)-1]; }
		inline unsigned int size(const int order) const { return Diag[order/2].size(); }
		inline T& FindH(const int order, const int m, const int n) { return Diag[order/2 + m*n/2- 1][Access::PosFromMN(m,n)-1]; }
		inline T& FindH(const int order, const int loc) { return FindH(order, Access::mAtLoc(loc), Access::nAtLoc(loc)); }
		inline int LengthOfThisOrder(const int order) { return rowLengths[order/2]; }

		inline static void StartThread(Hmn_c<T>* Hmn, const int startingPos, const int endingPos, const int order, T& temp1){
			Hmn->ThreadFillHmn(startingPos, endingPos, order, temp1);
		}
		
		static void FillHmn(Hmn_c<T>& Hmn);
};

template<class T>
Hmn_c<T>::Hmn_c(Cpqmn_c<T>* Cpqmn, const std::vector<T> hp, const int maxOrder, const bool complexArgs): Cpqmn(Cpqmn), hp(hp), HIsComplex(complexArgs)
{
/*	if(GetAvailableMemory() < (long)maxOrder*numberOfMN*precision/16){
		swapping = true;
		std::cout << "I think the available memory, " << GetAvailableMemory() << " isn't enough to store our Hmn at " << (long)maxOrder*numberOfMN*precision/16 << "." << std::endl;
	} else {
		swapping = false;
		std::cout << "I think the available memory, " << GetAvailableMemory() << " is enough to store our Hmn at " << (long)maxOrder*numberOfMN*precision/16 << "." << std::endl;
	}*/
	static_assert(std::is_same<T,mpfr::mpreal>::value || std::is_same<T,std::complex<mpfr::mpreal>>::value, "Hmn_c must be instantiated as Hmn_c<mpreal> or Hmn_c<complex<mpreal>>.");
	Diag.resize(maxOrder/2);
	rowLengths.resize(maxOrder/2);
	int diagSize = 0;
	for(int order = 0; order <= maxOrder-2; order+=2){
		diagSize += Access::MultOfMN(order+2);
		Diag[order/2].resize(diagSize);
		rowLengths[(maxOrder-order)/2-1] = diagSize;
	}
	for(int pos = 1; pos <= diagSize; ++pos) FindH(0, pos-1) = 1;
	if(HIsComplex){
		complexH.resize(hp.size());
		for(unsigned int i = 0; i < complexH.size(); ++i){
			complexH[i].resize(maxOrder/2+1);
			complexH[i][0] = 1;
		}
	} else {
		realH.resize(hp.size());
		for(unsigned int i = 0; i < realH.size(); ++i){
			realH[i].resize(maxOrder/2+1);
			realH[i][0] = 1;
		}
	}
//	std::cout << "After allocating, there is now " << GetAvailableMemory() << " available." << std::endl;
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
void Hmn_c<T>::FillHmn(Hmn_c<T>& Hmn){
	const int maxOrder = Hmn.Diag.size()*2;
	std::vector<std::thread> thread;
	thread.resize(maxThreads);
	std::vector<T> temp1;
	temp1.resize(maxThreads);
	int numThreads, posThisOrder, posPerThread, extraPos;
	float progress, totalComputations=0.0f;
	if(showProgressBar){
		progress = 0.0f;
		for(int order = 2; order < maxOrder; order+=2) totalComputations += Hmn.Diag[(maxOrder-order)/2-1].size()*Hmn.Diag[order/2-1].size();
		std::cout << "\r";
		DrawProgressBar(progress);
	}
	for(int order = 2; order < maxOrder; order+=2){	
		posThisOrder = Hmn.LengthOfThisOrder(order);
		numThreads = std::min(maxThreads, posThisOrder);
		posPerThread = posThisOrder / numThreads;
		extraPos = posThisOrder % numThreads;
		thread[0] = std::thread(Hmn_c::StartThread, &Hmn, 1, posPerThread + extraPos, order, std::ref(temp1[0]));
		for(int i=2; i<=numThreads; ++i){	
			thread[i-1] = std::thread(Hmn_c::StartThread, &Hmn, (i-1)*posPerThread + extraPos + 1, i*posPerThread + extraPos, order, std::ref(temp1[i-1]));
		}
		for(int i=1; i<= numThreads; ++i){
			thread[i-1].join();
		}
		FillH(Hmn, order, temp1[0]);
		if(showProgressBar){
			progress += posThisOrder*Hmn.Diag[order/2-1].size()/totalComputations;
			DrawProgressBar(progress);
		}
		Hmn.Diag[order/2-1].clear();
		Hmn.Cpqmn->DoneWithOrder(order);
	}
	FillH(Hmn, maxOrder, temp1[0]);
	Hmn.Diag[maxOrder/2-1].clear();
	if(showProgressBar){
		progress = 1.0f;
		DrawProgressBar(progress);
		std::cout << std::endl;
	}
	return;
}

template<class T>
void Hmn_c<T>::ThreadFillHmn(const int startingPos, const int endingPos,
		const int order, T& temp1){
	for(int pos = startingPos; pos <= endingPos; ++pos){
		for(unsigned int scanPos = 1; scanPos <= Diag[order/2-1].size(); ++scanPos){
			temp1 = (*Cpqmn)[pos-1][scanPos-1]*Diag[order/2-1][scanPos-1];
			FindH(order, Access::mAtLoc(pos-1), Access::nAtLoc(pos-1)) += temp1;
		}
	}
	return;
}

template<>
inline void Hmn_c<mpfr::mpreal>::FillH(
		Hmn_c<mpfr::mpreal>& Hmn, const int order, 
		mpfr::mpreal temp){
	for(unsigned int i = 0; i < Hmn.realH.size(); ++i){
		for(unsigned int scanPos = 1; scanPos <= Hmn.size(order-2); ++scanPos){
			temp = Hmn.Cpqmn->CFromHp(Hmn.hp[i], scanPos)*Hmn.at(order-2, scanPos-1);
			Hmn.realH[i][order/2] += temp;
		}
	}
	return;
}

template<>
inline void Hmn_c<std::complex<mpfr::mpreal>>::FillH(
		Hmn_c<std::complex<mpfr::mpreal>>& Hmn, const int order, 
		std::complex<mpfr::mpreal> temp){
	if(Hmn.HIsComplex){
		for(unsigned int i = 0; i < Hmn.complexH.size(); ++i){
			for(unsigned int scanPos = 1; scanPos <= Hmn.size(order-2); ++scanPos){
				temp = Hmn.Cpqmn->CFromHp(Hmn.hp[i], scanPos)*Hmn.at(order-2, scanPos-1);
				Hmn.complexH[i][order/2] += temp;
			}
		}
	} else {
		for(unsigned int i = 0; i < Hmn.realH.size(); ++i){
			for(unsigned int scanPos = 1; scanPos <= Hmn.size(order-2); ++scanPos){
				temp = Hmn.Cpqmn->CFromHp(Hmn.hp[i], scanPos)*Hmn.at(order-2, scanPos-1);
				Hmn.realH[i][order/2] += temp.real();
			}
		}
	}
	return;
}
} // namespace virasoro
#endif
