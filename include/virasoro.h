#ifndef VIRASORO_H_
#define VIRASORO_H_

#include <chrono>		// timers
#include <iostream>		// cout
#include <fstream>		// file output
#include <string>		// std::string
#include <vector>		// std::vector

#ifdef HAVE_WSTP
#include "wstp.h"
#endif
#include "mpreal.h"
#include "mpcomplex.h"
#include "runfile.h"
#include "access.h"
#include "virasoro.h"
#include "cpqmn.h"
#include "hmn.h"

namespace virasoro {

typedef std::chrono::high_resolution_clock Clock;

int core(int argc, char** argv, const bool wolfram);
std::vector<std::string> CollectArgs(int argc, char** argv);
void ReadDefaults(const std::string filename, const bool quiet);
void CreateConfigFile(const std::string filename);
std::string ParseOptions(std::vector<std::string>& args);
void DoOptions(const std::string& options, const bool quiet);

std::string NameOutputFile(const Runfile_c& runfile);
bool ParamsReal(const std::vector<std::complex<mpfr::mpreal>>& runVec);
void CheckRealityAndRun(const std::vector<std::complex<mpfr::mpreal>>& runVec, const int maxOrder, const std::string outputName, const int bGiven);
int DoRuns(const Runfile_c& runfile, const std::string options);
void ShowTime(const std::string computationName, const std::chrono::time_point<std::chrono::high_resolution_clock> timeStart);

std::string to_string(const mpfr::mpreal& N, int digits);
std::string to_string(const std::complex<mpfr::mpreal>& N, int digits, int base = 10);

template<class T>
void ConvertInputs(T& bsq, T& invBsq, T& llsq, T& lhsq, const T& c, const T& hl, const T& hh, T& temp1, T& temp2){
	temp1 = c*c;
	temp2 = c*26;
	temp1 -= temp2;
	temp1 += 25;
	temp1 = mpfr::sqrt(temp1);
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

/* !! this should actually be a static Hmn_c function that fills one order of H
template<class T>
void FillH(T* H, const Hmn_c<T>* Hmn, const Cpqmn_c<T>* Cpqmn, const T hp, const unsigned short int maxOrder){
	T temp1;
	for(int order = 2; order <= maxOrder; order+=2){
		for(unsigned int scanPos = 1; scanPos <= Hmn->size(order-2); ++scanPos){
			temp1 = hp - Cpqmn->hpmn[scanPos-1];
			temp1 = Cpqmn->Rmn[scanPos-1]/temp1; // Rmn includes the factor 16^mn
			temp1 *= Hmn->at(order-2, scanPos-1);
			H[order/2] += temp1;
		}
	}
	return;
}*/

template<class T, class U>
void DisplayH(const std::vector<U> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder){
	std::cout << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2]/*.real()*/, 10) << std::endl;
	}
}

// An empty outputName means a single run; a filled one is a multirun, which prints fewer words.
template<class T, class U>
void WriteH(const std::vector<U> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder, const std::string outputName){
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
#ifdef HAVE_WSTP
template<class T, class U>
void WSTPOut(const std::vector<U> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder){
	WSPutFunction(stdlink, "List", 5);
	WSPutString(stdlink, to_string(c,0).c_str());
	WSPutString(stdlink, to_string(hl,0).c_str());
	WSPutString(stdlink, to_string(hh,0).c_str());
	WSPutString(stdlink, to_string(hp,0).c_str());
	std::string moString(std::to_string(maxOrder));
	WSPutString(stdlink, moString.c_str()); 
	WSPutFunction(stdlink, "List", maxOrder/2 + 1);
	for(unsigned int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		WSPutString(stdlink, to_string(H[orderBy2],0).c_str());
	}
	return;
}
#endif

template<class T>
void FindCoefficients(std::vector<T> runVector, unsigned short int maxOrder, const std::string outputName, const int bGiven, const bool complexH){
#ifdef VERBOSE_DEBUG
	std::cout << "Beginning run with";
	std::cout << " c=" << to_string(runVector[0], 4);
	std::cout << " hl=" << to_string(runVector[1], 4);
	std::cout << " hh=" << to_string(runVector[2], 4);
	std::cout << " hp=";
	for(unsigned int i = 4; i <= runVector.size(); ++i) std::cout << to_string(runVector[i-1], 4) << ",";
	std::cout << "\b " << std::endl;
#endif
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	if(showProgressBar){
		std::cout << "Computing prefactors...";
		std::cout.flush();
	}
	T bsq, invBsq, llsq, lhsq, temp1, temp2;
	if(bGiven == 1){
		bsq = runVector[0];
		bsq *= runVector[0];
		invBsq = 1/bsq;
		runVector[0] = static_cast<std::complex<mpfr::mpreal>>(13) + 6*(bsq + invBsq);
	}
	if(bGiven == 2){
		bsq = runVector[0];
		invBsq = 1/bsq;
		runVector[0] = static_cast<std::complex<mpfr::mpreal>>(13) + 6*(bsq + invBsq);
	}
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);

	maxOrder -= maxOrder%2;
	Access::Populate(maxOrder);
	
	Cpqmn_c<T> Cpqmn(bsq, invBsq, maxOrder);
	
	Cpqmn_c<T>::FillHpmn(Cpqmn);
	Cpqmn_c<T>::FillRmn(Cpqmn, &llsq, &lhsq);
	if(outputName == "__WSTP" || outputName == "__MATHEMATICA"){
		maxOrder = Cpqmn_c<T>::CheckForDivergences(Cpqmn, showProgressBar, true);
	} else {
		maxOrder = Cpqmn_c<T>::CheckForDivergences(Cpqmn, showProgressBar, false);
	}
	if(maxOrder <= 2) return;
	Cpqmn_c<T>::FillCpqmn(Cpqmn);

	// combine Rmn and hpmn into computation of H
	std::vector<T> hpvec;
	for(unsigned int i = 4; i <= runVector.size(); ++i) hpvec.push_back(runVector[i-1]);
	Hmn_c<T> Hmn(&Cpqmn, hpvec, maxOrder, complexH);
	Hmn_c<T>::FillHmn(Hmn);
	
	// record and display coefficients H by q order
	if(complexH){
		for(unsigned int i = 4; i <= runVector.size(); ++i){
			if(outputName.empty() || outputName == "__CONSOLE") DisplayH(Hmn.complexH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
			if(outputName != "__CONSOLE" && outputName != "__WSTP") WriteH(Hmn.complexH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder, outputName);
#ifdef HAVE_WSTP
			if(outputName == "__WSTP") WSTPOut(Hmn.complexH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
#endif
		}
	} else {
		for(unsigned int i = 4; i <= runVector.size(); ++i){
			if(outputName.empty() || outputName == "__CONSOLE") DisplayH(Hmn.realH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
			if(outputName != "__CONSOLE" && outputName != "__WSTP") WriteH(Hmn.realH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder, outputName);
#ifdef HAVE_WSTP
			if(outputName == "__WSTP") WSTPOut(Hmn.realH[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
#endif
		}
	}
	return;
}
} // namespace virasoro
#endif
