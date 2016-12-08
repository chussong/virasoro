#ifndef VIRASORO_H_
#define VIRASORO_H_
//#define _CRT_DISABLE_PERFCRIT_LOCKS			requires removing iostream and possibly extra static linking

#include <cstdlib>		// atoi
#include <chrono>		// timers
#include <iostream>		// cout
#include <string>		// std::string
#include <vector>		// std::vector
#include <fstream>		// file output
#include <stdio.h>		// fgetc
#include <gmpxx.h>		// mpf_class
#include <thread>		// std::thread
#include <tuple>		// std::tuple

#include "mpfc_class.h"
#include "cpqmn.h"
#include "hmn.h"
#include "runfile.h"

extern int maxThreads;
extern int precision;
extern mpf_class tolerance;
extern bool showProgressBar;

typedef std::chrono::high_resolution_clock Clock;

std::vector<std::string> CollectArgs(int argc, char** argv);

void ReadDefaults(std::string filename);

void CreateConfigFile(std::string filename);

std::string ParseOptions(std::vector<std::string> &args);

int RunFromFile(std::string filename, const std::string options);

int RunFromTerminal(std::vector<std::string> args, const std::string options);

int ExecuteRunfile(Runfile_c runfile, std::string options);

void DebugPrintRunVector(const mpf_class* runVector, const std::vector<mpf_class> hp, const unsigned short int maxOrder);

int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);

void ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart);

std::string to_string(const mpf_class N, int digits);

template<class T>
void FindCoefficients(std::vector<T> runVector, unsigned short int maxOrder, const std::string outputName, const int bGiven){
/*	std::cout << "Beginning run with";
	std::cout << " c=" << to_string(runVector[0], 4);
	std::cout << " hl=" << to_string(runVector[1], 4);
	std::cout << " hh=" << to_string(runVector[2], 4);
	std::cout << " hp=";
	for(unsigned int i = 4; i <= runVector.size(); ++i) std::cout << to_string(runVector[i-1], 4) << ",";
	std::cout << "\b " << std::endl;*/
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	if(showProgressBar){
		std::cout << "Computing prefactors...";
		std::cout.flush();
	}
	T bsq, invBsq, llsq, lhsq, temp1, temp2;
	if(bGiven == 1){
		bsq = runVector[0]*runVector[0];
		bsq *= runVector[0];
		invBsq = 1/bsq;
		runVector[0] = 13 + 6*(bsq + invBsq);
	}
	if(bGiven == 2){
		bsq = runVector[0];
		invBsq = 1/bsq;
		runVector[0] = 13 + 6*(bsq + invBsq);
	}
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);

	int* mnLocation = new int[maxOrder]; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	int* mnMultiplicity = new int[maxOrder]();	/* number of mn combinations giving i+1 = m*n */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int* mTable = new unsigned short int[numberOfMN]();
	unsigned short int* nTable = new unsigned short int[numberOfMN]();
	int* mnLookup = new int[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, mnLocation, mnMultiplicity, maxOrder);
	
	Cpqmn_c<T> Cpqmn(&bsq, &invBsq, numberOfMN, maxOrder, mTable, nTable, mnLookup);
	Cpqmn.FillHpmn();
	
	auto time1 = Clock::now();
	Cpqmn.FillRmn(&llsq, &lhsq);
	auto time2 = Clock::now();
	CheckForDivergences(&Cpqmn, maxOrder, mnLocation, mnMultiplicity);
	if(maxOrder <= 2) return;
	Cpqmn.FillCpqmn();

	// combine Rmn and hpmn into computation of H
	Hmn_c<T> Hmn(&Cpqmn, numberOfMN, maxOrder, mnLocation, mnMultiplicity, mnLookup);
	time1 = Clock::now();
	Hmn.FillHmn();
	time2 = Clock::now();
	
	// corral H by q order and display coefficients
	T* H = new T[maxOrder/2+1];
	for(unsigned int i = 4; i <= runVector.size(); ++i){
		H[0] = 1;
		FillH(H, &Hmn, &Cpqmn, runVector[i-1], mnLocation, mnMultiplicity, maxOrder);
		if(outputName.empty() || outputName == "__CONSOLE") DisplayH(H, runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
		if(outputName != "__CONSOLE") WriteH(H, runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder, outputName);
	}
	delete[] mnLocation;
	delete[] mnMultiplicity;
	delete[] mTable;
	delete[] nTable;
	delete[] mnLookup;
	delete[] H;
	return;
}

template<class T>
void CheckForDivergences(const Cpqmn_c<T>* Cpqmn, unsigned short int &maxOrder, const int* mnLocation, const int* mnMultiplicity){
	int oldMax = maxOrder;
	for(int pq = 2; pq <= maxOrder-2; pq+=2){
		for(int pos = mnLocation[pq-1]; pos < mnLocation[pq-1] + mnMultiplicity[pq-1]; ++pos){
			for(int mn = pq; mn <= maxOrder-2; mn+=2){
				for(int scanPos = mnLocation[mn-1]; scanPos < mnLocation[mn-1] + mnMultiplicity[mn-1]; ++scanPos){
					if(abs(Cpqmn->hpmn[pos-1] + pq - Cpqmn->hpmn[scanPos-1]) < tolerance && maxOrder > std::max(mn-2,pq-2)) maxOrder = std::max(mn-2, pq-2);
				}
			}
			if(Cpqmn->Amn[pos-1] == 0 && maxOrder > pq-2){
				maxOrder = pq-2;
				if(maxOrder > 2) printf("Stopping this run at order %i because Amn diverges above this.\n", maxOrder);
				if(maxOrder <= 2) printf("Skipping this run because Amn diverges immediately.\n");
				return;
			}
		}
 	}
	maxOrder = maxOrder - (maxOrder%2);
	if(maxOrder < oldMax){
		if(maxOrder > 2) printf("Stopping this run at order %i because the coefficients diverge above this.\n",maxOrder);
		if(maxOrder <= 2) printf("Skipping this run because the coefficients diverge immediately.\n");
	}
	return;
}

template<class T>
void ConvertInputs(T& bsq, T& invBsq, T& llsq, T& lhsq, const T& c, const T& hl, const T& hh, T& temp1, T& temp2){
	temp1 = c*c;
	temp2 = c*26;
	temp1 -= temp2;
	temp1 += 25;
	temp1 = sqrt(temp1);
	temp1 = c - temp1;
	temp1 -= 13;
	bsq = temp1/12;
	invBsq = 1/bsq;

	temp1 = c - 1;
	temp1 /= 24;
	llsq = hl - temp1;
	temp1 = c - 1;
	temp1 /= 24;
	lhsq = hh - temp1;
}

template<class T>
void FillH(T* H, const Hmn_c<T>* Hmn, const Cpqmn_c<T>* Cpqmn, const T hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	T temp1;
	for(int order = 2; order <= maxOrder; order+=2){
		H[order/2] = 0;
		for(int power = 2; power <= order; power+=2){
			for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
				temp1 = hp - Cpqmn->hpmn[scanPos-1];
				temp1 = Cpqmn->Rmn[scanPos-1]/temp1;
				temp1 <<= 4*power;		// fast multiplication by 2^(4*power)
//				std::cout << "H[" << order << "] += (" << to_string(temp1, 4) << ")*(";
				temp1 *= Hmn->Hmn[(order-power)/2][scanPos-1];
//				std::cout << to_string(Hmn->Hmn[(order-power)/2][scanPos-1], 4) << ") = (" << to_string(temp1, 4) << ")" << std::endl;
				H[order/2] += temp1;
			}
		}
	}
	return;
}

template<class T>
void DisplayH(const T* H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder){
	std::cout << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
	}
}

// An empty outputName means a single run; a filled one is a multirun, which prints fewer words.
template<class T>
void WriteH(const T* H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder, const std::string outputName){
	std::ofstream outputFile;
	std::string filename = outputName;
	if(outputName == "__MATHEMATICA"){
		std::cout << "{" << to_string(c, 0) << "," << to_string(hl, 0) << "," << to_string(hh, 0) << "," << to_string(hp, 0) << "," << maxOrder << "},{1";
		for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
			std::cout << "," << to_string(H[orderBy2], 0);
		}
		std::cout << "}";
		return;
	}
	outputFile.open (filename, std::ios_base::app | std::ios_base::out);
	if(outputName.empty()){
		outputFile << "Given the parameters" << std::endl;
		outputFile << "c = " << to_string(c, 0) << ", h_L = " << to_string(hl, 0) << ", h_H = " << to_string(hh, 0) << ", h_p = " << to_string(hp, 0) << std::endl;
		outputFile << "the Virasoro block coefficients are as follows:" << std::endl;
		for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
			outputFile << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
		}
	} else {
		outputFile << "{" << to_string(c, 0) << "," << to_string(hl, 0) << "," << to_string(hh, 0) << "," << to_string(hp, 0) << "," << maxOrder << "}" << std::endl;
	}
	outputFile << "{1";
	for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
		outputFile << "," << to_string(H[orderBy2], 0);
	}
	outputFile << "}" << std::endl;
	outputFile.close();
}

#endif
