#include "hmn.h"

Hmn_t::Hmn_t(Cpqmn_t* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup): numberOfMN(numberOfMN), maxOrder(maxOrder), mnLocation(mnLocation), mnMultiplicity(mnMultiplicity), mnLookup(mnLookup), Cpqmn(Cpqmn)
{
	Hmn = new mpf_class*[maxOrder/2];
	Hmn[0] = new mpf_class[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos) Hmn[0][pos-1] = 1;
	mnAtThisOrder = new int[maxOrder/2];
	for(int order = 2; order <= maxOrder-2; order+=2){
		mnAtThisOrder[order/2-1] = 0;
		for(int mn = 2; mn <= maxOrder-order; mn+=2){
			mnAtThisOrder[order/2-1] += mnMultiplicity[mn-1];
		}
		Hmn[order/2] = new mpf_class[mnAtThisOrder[order/2-1]];
	}
}

Hmn_t::~Hmn_t(){
	delete[] Hmn[0];
	for(int order = 2; order <= maxOrder-2; order +=2){
		delete[] Hmn[order/2];
	}
	delete[] Hmn;
	delete[] mnAtThisOrder;
}

void Hmn_t::FillHmn(){	
	std::thread* thread = new std::thread[maxThreads];
//	mpf_class temp1[maxThreads], temp2[maxThreads], hpTemp[maxThreads];
	mpf_class* temp1 = new mpf_class[maxThreads];
	mpf_class* temp2 = new mpf_class[maxThreads];
	mpf_class* hpTemp = new mpf_class[maxThreads];
	int numThreads;
	for(int order = 2; order < maxOrder; order+=2){	
		numThreads = std::min(maxThreads,(maxOrder-order)/2);
		thread[0] = std::thread(Hmn_t::StartThread, this, 2, (maxOrder-order)/numThreads + (maxOrder-order)%numThreads, order, std::ref(temp1[0]), std::ref(temp2[0]), std::ref(hpTemp[0]));
		for(int i=2; i<=numThreads; ++i){	
			thread[i-1] = std::thread(Hmn_t::StartThread, this, (i-1)*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads + 1, i*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads, order, std::ref(temp1[i-1]), std::ref(temp2[i-1]), std::ref(hpTemp[i-1]));
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

void Hmn_t::ThreadFillHmn(const int startingMN, const int endingMN, const int order, mpf_class& temp1, mpf_class& temp2, mpf_class& hpTemp){
	for(int mn = startingMN + startingMN%2; mn <= endingMN; mn+=2){
		for(int pos = mnLocation[mn-1]; pos <= mnLocation[mn-1] + mnMultiplicity[mn-1] - 1; ++pos){
			hpTemp = Cpqmn->hpmn[pos-1] + mn;
			for(int power = 2; power <= order; power+=2){
				for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
					temp1 = Cpqmn->Rmn[scanPos-1];
					temp1 *= Hmn[(order-power)/2][scanPos-1];
					temp1 *= powOverflow[power/256];
					temp2 = 16;
					mpf_pow_ui(temp2.get_mpf_t(), temp2.get_mpf_t(), power%256);
					temp1 *= temp2;
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
