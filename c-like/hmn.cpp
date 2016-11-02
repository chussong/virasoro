#include "hmn.h"

Hmn_t::Hmn_t(Cpqmn_t* Cpqmn, const int numberOfMN, const unsigned short int maxOrder, const int* mnLocation, const int* mnMultiplicity, const int* mnLookup): numberOfMN(numberOfMN), maxOrder(maxOrder), mnLocation(mnLocation), mnMultiplicity(mnMultiplicity), Cpqmn(Cpqmn)
{
	Hmn = new mpf_t*[maxOrder/2];
	Hmn[0] = new mpf_t[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos) mpf_init_set_ui(Hmn[0][pos-1], 1);
	mnAtThisOrder = new int[maxOrder/2];
	for(int order = 2; order <= maxOrder-2; order+=2){
		mnAtThisOrder[order/2-1] = 0;
		for(int mn = 2; mn <= maxOrder-order; mn+=2){
			mnAtThisOrder[order/2-1] += mnMultiplicity[mn-1];
		}
		Hmn[order/2] = new mpf_t[mnAtThisOrder[order/2-1]];
		for(int pos = 1; pos <= mnAtThisOrder[order/2-1]; ++pos) mpf_init(Hmn[order/2][pos-1]);
	}
}

Hmn_t::~Hmn_t(){
	for(int pos = 1; pos <= numberOfMN; ++pos) mpf_clear(Hmn[0][pos-1]);
	delete[] Hmn[0];
	for(int order = 2; order <= maxOrder-2; order +=2){
		for(int pos = 1; pos <= mnAtThisOrder[order/2-1]; ++pos){
			mpf_clear(Hmn[order/2][pos-1]);
		}
		delete[] Hmn[order/2];
	}
	delete[] Hmn;
	delete[] mnAtThisOrder;
	Hmn = nullptr;
	mnAtThisOrder = nullptr;
}

void Hmn_t::FillHmn(){
	std::thread thread[maxThreads];
	mpf_t temp1[maxThreads], temp2[maxThreads], hpTemp[maxThreads];
	for(int i=1; i <= maxThreads; ++i) mpf_inits(temp1[i-1], temp2[i-1], hpTemp[i-1], NULL);
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
	for(int i=1; i <= maxThreads; ++i) mpf_clears(temp1[i-1], temp2[i-1], hpTemp[i-1], NULL);
	return;
}

void Hmn_t::ThreadFillHmn(const int startingMN, const int endingMN, const int order, mpf_t& temp1, mpf_t& temp2, mpf_t& hpTemp){
	for(int mn = startingMN + startingMN%2; mn <= endingMN; mn+=2){
		for(int pos = mnLocation[mn-1]; pos <= mnLocation[mn-1] + mnMultiplicity[mn-1] - 1; ++pos){
			mpf_add_ui(hpTemp, Cpqmn->hpmn[pos-1], mn);
			for(int power = 2; power <= order; power+=2){
				for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
					mpf_sub(temp1, hpTemp, Cpqmn->hpmn[scanPos-1]);
					mpf_div(temp1, Cpqmn->Rmn[scanPos-1], temp1);
					mpf_mul(temp1, temp1, Hmn[(order-power)/2][scanPos-1]);
					mpf_mul(temp1, temp1, powOverflow[power/256]);
					mpf_set_ui(temp2, 16);
					mpf_pow_ui(temp2, temp2, power%256);
					mpf_mul(temp1, temp1, temp2);
					mpf_add(Hmn[order/2][pos-1], Hmn[order/2][pos-1], temp1);
				}
			}
		}
	}
	return;
}
