#ifndef RUNFILE_H_ 
#define RUNFILE_H_

#include <iostream>
#include <fstream>
#include <algorithm>	// std::swap
#include <string>
#include <vector>
#include <tuple>
#ifdef HAVE_WSTP_
#include "wstp.h"
#endif
#include "mpreal.h"
#include "mpcomplex.h"
#include "access.h"
#include "cpqmn.h"
#include "hmn.h"
//#include "config.h" !! is this actually necessary?

typedef std::chrono::high_resolution_clock Clock;

static int 				static_maxOrder;
static std::vector<int> static_mnLocation;
static std::vector<int> static_mnLookup;
static std::vector<int> static_mTable;
static std::vector<int> static_nTable;

//inline int mnpos(const int m, const int n) { std::cout << "Accessing element " << m*n-1 << " of mnpos, which is size " << static_mnLocation.size() << ".\n"; return static_mnLocation[m*n-1]; }
inline int mnpos(const int m, const int n) { return static_mnLookup[(m-1)*static_maxOrder + n-1]; }
inline int MAt(int loc)                    { return static_mTable[loc];       }
inline int NAt(int loc)                    { return static_nTable[loc];       }

std::string to_string(const mpfr::mpreal N, int digits);
std::string to_string(const std::complex<mpfr::mpreal> N, int digits, int base = 10);

class Runfile_c{
	std::string filename;
	std::vector<std::string> lines;

	int maxThreads = 8;				// Maximum number of simultaneous threads
	int precision = 768;			// Precision of mpfr::mpreal and std::complex<mpfr::mpreal> in bits
	mpfr::mpreal tolerance = 1e-10;		// Smaller than this is taken to be 0 for comparisons
	bool showProgressBar = true;	// Show progress bar during FillHmn()

	public:
		std::vector<std::vector<std::complex<mpfr::mpreal>>> runs;
		std::vector<int> maxOrders;

		Runfile_c();
		Runfile_c(const Runfile_c& other);
		Runfile_c(Runfile_c&& other);
		Runfile_c(const char* filename);
		Runfile_c(const std::string filename);
		Runfile_c(const std::vector<std::string> line);
		~Runfile_c();
		void swap(Runfile_c& first, Runfile_c& second);

		Runfile_c& operator=(Runfile_c&& v);
		Runfile_c& operator=(Runfile_c v);
		Runfile_c& operator=(const char* &filename);
		Runfile_c& operator=(const std::string& filename);
		Runfile_c& operator=(const std::vector<std::string> line);
		
		void SetMaxThreads(int newMax);
		void SetPrecision(int newPrec);
		void SetTolerance(mpfr::mpreal newTolerance);
		void SetProgressBar(bool newProgressBar);
		int NumberOfRuns();

		int ReadRunfile();
		int Execute(std::string options);

		int Expand();
		int ExpandBraces(const int param);
		std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> ParseBraces(std::string insideBraces);
		std::tuple<size_t, size_t> FindNthParameter(const std::string line, const int param);

		int ExpandRelativeEqns(const int param);
		std::tuple<std::complex<mpfr::mpreal>, int> ParseRelativeEqn(std::string equation, std::string relTo);
		std::complex<mpfr::mpreal> RelativeMPF(std::string firstHalf, std::string equation);
		std::string FindBaseNumber(std::string sourceString, const int paramNumber);

		int RunCompare(std::vector<std::complex<mpfr::mpreal>> run1, std::vector<std::complex<mpfr::mpreal>> run2);

		std::string NameOutputFile();

		int EnumerateMN (int* mnLocation, int* mnMultiplicity,  unsigned short int maxOrder);
		void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder);
		void ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart);

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
			int* mnLocation = new int[maxOrder]; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
			int* mnMultiplicity = new int[maxOrder]();	/* number of mn combinations giving i+1 = m*n */
			int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
			unsigned short int* mTable = new unsigned short int[numberOfMN]();
			unsigned short int* nTable = new unsigned short int[numberOfMN]();
			int* mnLookup = new int[maxOrder*maxOrder];
			FillMNTable(mnLookup, mTable, nTable, mnLocation, mnMultiplicity, maxOrder);
			Access::Populate(maxOrder);
			
			Cpqmn_c<T> Cpqmn(bsq, invBsq, maxOrder);
			Cpqmn.FillHpmn();
			
			auto time1 = Clock::now();
			Cpqmn.FillRmn(&llsq, &lhsq);
			auto time2 = Clock::now();
			if(outputName == "__WSTP" || outputName == "__MATHEMATICA"){
				maxOrder = Cpqmn_c<T>::CheckForDivergences(Cpqmn, showProgressBar, true);
			} else {
				maxOrder = Cpqmn_c<T>::CheckForDivergences(Cpqmn, showProgressBar, false);
			}
			if(maxOrder <= 2) return;
			Cpqmn.FillCpqmn();

			// combine Rmn and hpmn into computation of H
			std::vector<T> hpvec;
			for(unsigned int i = 4; i <= runVector.size(); ++i) hpvec.push_back(runVector[i-1]);
			Hmn_c<T> Hmn(&Cpqmn, hpvec, maxOrder);
			time1 = Clock::now();
			Hmn_c<T>::Fill(Hmn);
			time2 = Clock::now();
			
			// corral H by q order and display coefficients
//			T* H = new T[maxOrder/2+1];
			for(unsigned int i = 4; i <= runVector.size(); ++i){
/*				H[0] = 1;
				FillH(H, &Hmn, &Cpqmn, runVector[i-1], maxOrder);*/
				if(outputName.empty() || outputName == "__CONSOLE") DisplayH(Hmn.H[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
				if(outputName != "__CONSOLE" && outputName != "__WSTP") WriteH(Hmn.H[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder, outputName);
#ifdef HAVE_WSTP_
				if(outputName == "__WSTP") WSTPOut(Hmn.H[i-4], runVector[0], runVector[1], runVector[2], runVector[i-1], maxOrder);
#endif
			}
			delete[] mnLocation;
			delete[] mnMultiplicity;
			delete[] mTable;
			delete[] nTable;
			delete[] mnLookup;
//			delete[] H;
			return;
		}

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

		// ~~ this should actually be a static Hmn_c function that fills one order of H
		template<class T>
		void FillH(T* H, /*const*/ Hmn_c<T>* Hmn, const Cpqmn_c<T>* Cpqmn, const T hp, const unsigned short int maxOrder){
			T temp1;
			for(int order = 2; order <= maxOrder; order+=2){
				for(unsigned int scanPos = 1; scanPos <= Hmn->size(order-2); ++scanPos){
					temp1 = hp - Cpqmn->hpmn[scanPos-1];
					temp1 = Cpqmn->Rmn[scanPos-1]/temp1; // Rmn includes the factor 16^mn
					temp1 *= Hmn->at(order-2, scanPos-1);
					H[order/2] += temp1;
				}
			}
/*			for(int order = 2; order <= maxOrder; order+=2){
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
			}*/
			return;
		}

		template<class T>
		void DisplayH(const std::vector<T> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder){
			std::cout << "Given the parameters" << std::endl;
			std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
			std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
			for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
				std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
			}
		}

		// An empty outputName means a single run; a filled one is a multirun, which prints fewer words.
		template<class T>
		void WriteH(const std::vector<T> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder, const std::string outputName){
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
#ifdef HAVE_WSTP_
		template<class T>
		void WSTPOut(const std::vector<T> H, const T c, const T hl, const T hh, const T hp, const unsigned short int maxOrder){
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
};

#endif
