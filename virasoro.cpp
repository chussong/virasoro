#include "virasoro.h"

typedef std::chrono::high_resolution_clock Clock;

int maxThreads = 8;			// Maximum number of simultaneous threads
int precision = 512;				// Precision of mpf_class in bits
const mpf_class tolerance(1e-20);	// Smaller than this is taken to be 0 for comparisons
mpf_class* powOverflow;

int main(int argc, char** argv){
	auto programStart = Clock::now();
	std::string options = ParseOptions(argc, argv);
	mpf_set_default_prec(precision);
	int exitCode;
	switch(argc){
		default:	printf("Error: input either the name of a runfile or the five parameters c, hl, hh, hp, maxOrder\n");
					return EXIT_FAILURE;
		case 2:		if(options.find("m", 0) != std::string::npos){
						exitCode = RunFromFile(argv[1], options);
					} else {
						exitCode = RunFromFile(argv[1], options);
						ShowTime("Entire computation", programStart);
					}
					break;
		case 6:		if(options.find("m", 0) != std::string::npos){
						exitCode = RunFromTerminal(argv, options);
					} else {
						exitCode = RunFromTerminal(argv, options);
						ShowTime("Entire computation", programStart);
					}
					break;
	}

	return exitCode;
}

// some day I'll get around to turning argv into a std::vector<std::string> to begin with
std::string ParseOptions(int &argc, char** &argv){
	std::string options = "";
	char** newArgv;
	int newArgc = argc;
	bool* realArg = new bool[argc]();
	for(int i = 1; i < argc; ++i){
		if(argv[i][0] != '-'){
			realArg[i] = true;
			continue;
		} else if(strcmp(argv[i], "-m") == 0){
			options.append("m");
			--newArgc;
			realArg[i] = false;
		} else if(strcmp(argv[i], "-c") == 0){
			options.append("c");
			--newArgc;
			realArg[i] = false;
		} else if(strcmp(argv[i], "-b") == 0){
			options.append("b");
			--newArgc;
			realArg[i] = false;
		} else if(strcmp(argv[i], "-bb") == 0){
			options.append("bb");
			--newArgc;
			realArg[i] = false;
		} else if(strncmp(argv[i], "-p", 2) == 0){
			std::string precString = argv[i];
			precision = std::stoi(precString.substr(2), nullptr, 10);
			--newArgc;
			realArg[i] = false;
		} else if(strncmp(argv[i], "-t", 2) == 0){
			std::string threadString = argv[i];
			maxThreads = std::stoi(threadString.substr(2), nullptr, 10);
			--newArgc;
			realArg[i] = false;
		} else {
			realArg[i] = true;
		}
	}
	newArgv = new char*[newArgc];
	newArgv[0] = argv[0];
	int pos = 1;
	for(int i = 1; i < argc; ++i){
		if(realArg[i]){
			newArgv[pos] = argv[i];
			++pos;
		}
	}
	argc = newArgc;
	argv = newArgv;
	delete[] newArgv;
	delete[] realArg;
	return options;
}

int RunFromFile(char* filename, const std::string options){
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;	
	int bGiven = 0;
	if(options.find("b", 0) != std::string::npos) bGiven = 1;
	if(options.find("bb", 0) != std::string::npos) bGiven = 2;
	Runfile_c runfile(filename);
	if(runfile.ReadRunfile() <= 0) return -1;
	if(!wolframOutput){
		for(unsigned int i = 1; i <= runfile.runs.size(); ++i){
			std::cout << "Run " << i << ": ";
			for(int j = 1; j <= 3; ++j){
				std::cout << to_string(runfile.runs[i-1][j-1], 0) << " ";
			}
			if(runfile.runs[i-1].size() > 4) std::cout << "{";
			for(unsigned int j = 4; j <= runfile.runs[i-1].size(); ++j){
				std::cout << to_string(runfile.runs[i-1][j-1], 0) << ",";
			}
			if(runfile.runs[i-1].size() > 4){
				std::cout << "\b} ";
			} else {
				std::cout << "\b ";
			}
			std::cout << runfile.maxOrders[i-1];
			std::cout << std::endl;
		}
	}
	int highestMax = 0;
	for(unsigned int i = 1; i <= runfile.runs.size(); ++i){
		if(runfile.maxOrders[i-1] > highestMax) highestMax = runfile.maxOrders[i-1];
	}
	SetPowOverflow(highestMax);
	auto runStart = Clock::now();
	std::string outputName;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
		for(unsigned int run = 1; run < runfile.runs.size(); ++run){
			FindCoefficients(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			std::cout << ",";
		}
		FindCoefficients(runfile.runs[runfile.runs.size()-1], runfile.maxOrders[runfile.runs.size()-1], outputName, bGiven);
		std::cout << "}";
	} else {
		if(consoleOutput){
			outputName = "__CONSOLE";
		} else {
			outputName = NameOutputFile(filename);
			std::remove(outputName.c_str());
		}
		for(unsigned int run = 1; run <= runfile.runs.size(); ++run){
			runStart = Clock::now();
//			DebugPrintRunVector(runs[run-1], hp[run-1], maxOrders[run-1]);
			std::cout << "Beginning run " << run << " of " << runfile.runs.size() << "." << std::endl;
			FindCoefficients(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);
		}
	}
	return 0;
}

int RunFromTerminal(char** argv, const std::string options){
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;
	int bGiven = 0;
	if(options.find("b", 0) != std::string::npos) bGiven = 1;
	if(options.find("bb", 0) != std::string::npos) bGiven = 2;
	std::vector<mpf_class> runVector;
	for(int i = 1; i <= 4; ++i){
		std::cout << " \b";				// mystifyingly it segfaults without this
		runVector.emplace_back(argv[i]);
	}
	unsigned short int maxOrder = std::atoi(argv[5]);
	maxOrder -= (maxOrder % 2);
	SetPowOverflow(maxOrder);
//	DebugPrintRunVector(runVector, hp, maxOrder);
	std::string outputName;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
	} else if(consoleOutput) {
		outputName = "__CONSOLE";
	} else {
		outputName = NameOutputFile(nullptr);
		std::remove(outputName.c_str());
	}
	FindCoefficients(runVector, maxOrder, outputName, bGiven);
	if(wolframOutput) std::cout << "}";
	return 0;
}

void SetPowOverflow(unsigned short int maxOrder){
	powOverflow = new mpf_class[maxOrder/256+1];
	powOverflow[0] = 1;
	for(int i = 1; i <= maxOrder/256; ++i){
		powOverflow[i] = 16;
		mpf_pow_ui(powOverflow[i].get_mpf_t(), powOverflow[i].get_mpf_t(), 256*i);
	}
}

void DebugPrintRunVector(const mpf_class* runVector, const std::vector<mpf_class> hp, const unsigned short int maxOrder){
	for(unsigned int i = 1; i <= hp.size(); ++i){
		std::cout << "Pretend this is a run with c = " << to_string(runVector[0], 0) << ", hl = " << to_string(runVector[1], 0) << ", hh = " << to_string(runVector[2], 0) << ", hp = " << to_string(hp[i-1], 0) << ", maxOrder = " << maxOrder << "." << std::endl;
	}
}

void FindCoefficients(std::vector<mpf_class> runVector, unsigned short int maxOrder, const std::string outputName, const int bGiven){
	if(runVector[0] > 1 && runVector[0] < 25 && bGiven == 0){
		std::cout << "This run appears to have a c value of " << runVector[0] << ", which is between 1 and 25. This will result in a complex b^2 and currently can not be handled. If this is supposed to be a value of b or b^2 instead of c, run again with -b or -bb." << std::endl;
		return;
	}
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	mpf_class bsq, invBsq, llsq, lhsq, temp1, temp2;
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);
	if(bGiven == 1){
		bsq = runVector[0]*runVector[0];
		invBsq = 1/bsq;
		runVector[0] = 13 + 6*(bsq + invBsq);
	}
	if(bGiven == 2){
		bsq = runVector[0];
		invBsq = 1/bsq;
		runVector[0] = 13 + 6*(bsq + invBsq);
	}
	CheckForDivergences(&invBsq, maxOrder);
	if(maxOrder <= 2) return;

	int* mnLocation = new int[maxOrder]; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	int* mnMultiplicity = new int[maxOrder]();	/* number of mn combinations giving i+1 = m*n */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int* mTable = new unsigned short int[numberOfMN]();
	unsigned short int* nTable = new unsigned short int[numberOfMN]();
	int* mnLookup = new int[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, mnLocation, mnMultiplicity, maxOrder);
	
	Cpqmn_t Cpqmn(&bsq, &invBsq, numberOfMN, maxOrder, mTable, nTable, mnLookup);
	Cpqmn.FillHpmn();
	
	auto time1 = Clock::now();
	Cpqmn.FillRmn(&llsq, &lhsq);
	auto time2 = Clock::now();
	
	// combine Rmn and hpmn into computation of H
	Hmn_t Hmn(&Cpqmn, numberOfMN, maxOrder, mnLocation, mnMultiplicity, mnLookup);
	time1 = Clock::now();
	Hmn.FillHmn();
	time2 = Clock::now();
	
	// corral H by q order and display coefficients
	mpf_class* H = new mpf_class[maxOrder/2+1];
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

void CheckForDivergences(const mpf_class* bsq, unsigned short int &maxOrder){
	mpf_class temp;
	int lowest = maxOrder;
	for(int m = 1; m <= maxOrder-2; ++m){
		for(int n = 1+(m%2); m*n <= maxOrder-2; n+=1+(m%2)){
			for(int p = 1; p <= maxOrder-m*n; ++p){
				for(int q = 1+(p%2); p*q <= maxOrder-m*n; q+=1+(p%2)){
					temp = (*bsq)*(n-q) - (m+p);
					if(abs(temp) < tolerance){
						if(m*n < lowest) lowest = m*n;
//						if(p*q < lowest) lowest = p*q;
						continue;
					}
					temp = (*bsq)*(n+q) - (m-p);
					if(abs(temp) < tolerance){
						if(m*n < lowest) lowest = m*n;
//						if(p*q < lowest) lowest = p*q;
						continue;
					}
				}
			}
		}
	}
	if(*bsq == 2) lowest = 2; // this is a shitty workaround until I can fix the 1/0 in Amn
	lowest = lowest - (lowest%2);
	if(lowest < maxOrder){
		maxOrder = lowest;
		if(lowest > 2) printf("Stopping this run at order %i because the coefficients diverge above this.\n",lowest);
		if(lowest <= 2) printf("Skipping this run because the coefficients diverge immediately.\n");
	}
	return;
}

int EnumerateMN (int* mnLocation, int* mnMultiplicity, const unsigned short int maxOrder){
	int numberOfMN = 0;
	for(int m=1; m <= maxOrder; m+=2){
		for(int n=2; m*n <= maxOrder; n+=2){		// odd m, even n
			++mnMultiplicity[m*n-1];
			++numberOfMN;
		}
		for(int n=1; (m+1)*n <= maxOrder; ++n){		// even m
			++mnMultiplicity[(m+1)*n-1];
			++numberOfMN;
		}
	}
	mnLocation[1] = 1;
	for(int i = 4; i <= maxOrder; i+=2){
		mnLocation[i-1] = mnLocation[i-3] + mnMultiplicity[i-3];
	}
	return numberOfMN;
}

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	int pos;
	for(int m = 1; m <= maxOrder; m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){		// odd m, even n
			for(pos = mnLocation[m*n-1]; pos <= mnLocation[m*n-1]+mnMultiplicity[m*n-1]; ++pos){
				if(mTable[pos-1]==0){
					mTable[pos-1] = m;
					nTable[pos-1] = n;
					mnLookup[(m-1)*maxOrder + n-1] = pos;	
					break;
				}
			}
		}
		for(int n = 1; (m+1)*n <= maxOrder; ++n){	// even m, all n
			for(pos = mnLocation[(m+1)*n-1]; pos <= mnLocation[(m+1)*n-1]+mnMultiplicity[(m+1)*n-1]; ++pos){
				if(mTable[pos-1]==0){
					mTable[pos-1] = m+1;
					nTable[pos-1] = n;
					mnLookup[m*maxOrder + n-1] = pos;	
					break;
				}
			}
		}
	}
}

void ConvertInputs(mpf_class& bsq, mpf_class& invBsq, mpf_class& llsq, mpf_class& lhsq, const mpf_class& c, const mpf_class& hl, const mpf_class& hh, mpf_class& temp1, mpf_class& temp2){
	if(c <= 1 || c >= 25){
		temp1 = c*c;
		temp2 = c*26;
		temp1 -= temp2;
		temp1 += 25;
		temp1 = sqrt(temp1);
		temp1 = c - temp1;
		temp1 -= 13;
		bsq = temp1/12;
		invBsq = 1/bsq;
	}

	temp1 = c - 1;
	temp1 /= 24;
	llsq = hl - temp1;
	temp1 = c - 1;
	temp1 /= 24;
	lhsq = hh - temp1;
}

void FillH(mpf_class* H, const Hmn_t* Hmn, const Cpqmn_t* Cpqmn, const mpf_class hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	mpf_class temp1, temp2;
	for(int order = 2; order <= maxOrder; order+=2){
		for(int power = 2; power <= order; power+=2){
			for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
				temp1 = hp - Cpqmn->hpmn[scanPos-1];
				temp1 = Cpqmn->Rmn[scanPos-1]/temp1;
				temp1 *= Hmn->Hmn[(order-power)/2][scanPos-1];
				temp1 *= powOverflow[power/256];
				temp2 = 16;
				mpf_pow_ui(temp2.get_mpf_t(), temp2.get_mpf_t(), power%256);
				temp1 *= temp2;
				H[order/2] += temp1;
			}
		}
	}
	return;
}

void ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart){
	auto timeEnd = Clock::now();
	int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
	std::string unit = "ms";
	if(elapsed > 5000){
		elapsed = std::chrono::duration_cast<std::chrono::seconds>(timeEnd - timeStart).count();
		unit = "s";
		if(elapsed > 300){
			elapsed = std::chrono::duration_cast<std::chrono::minutes>(timeEnd - timeStart).count();
			unit = "m";
			if(elapsed > 300){
				elapsed = std::chrono::duration_cast<std::chrono::hours>(timeEnd - timeStart).count();
				unit = "hr";
			}
		}
	}
	std::cout << computationName << " took " << elapsed << unit << "." << std::endl;	
}

std::string to_string(const mpf_class N, int digits){
	if(digits < 0) digits = -digits;	
	mp_exp_t dotPos;
	std::string output = N.get_str(dotPos, 10, digits);
	if(output.empty()) output.append("0");
	if(digits > 0 && dotPos < digits){
		while(dotPos > (int)output.size()) output.append("0");
		if((N > 1 || N < 0) && dotPos < (int)output.size()) output.insert(abs(dotPos), ".");
	}
	if(digits == 0 && dotPos > 0){
		while(dotPos > (int)output.size()) output.append("0");
		if(dotPos < (int)output.size()){
			if(sgn(N) == -1) dotPos += 1;
			output.insert(dotPos, ".");
		}
	}
	if(digits > 0 && dotPos > digits){
		while((int)output.size() < digits) output.append("0");
		if(sgn(N) == 1) output.insert(1, ".");		
		if(sgn(N) == -1) output.insert(2, ".");
		output.append("*10^");
		output.append(std::to_string(dotPos));
	}
	if(N > 0 && N < 1) output.insert(0, "0.");
	return output;
}

std::string NameOutputFile(const char* runfileName){
	std::string filename;
	if(runfileName == nullptr){
		filename = "";
	} else {
		filename = runfileName;
		std::size_t delPos = filename.find(".txt");
		if(delPos != std::string::npos) filename.erase(delPos, 4);
		filename.append("_results.txt");
	}
	return filename;
}

void DisplayH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder){
	std::cout << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
	}
}

// An empty outputName means a single run; a filled one is a multirun, which prints fewer words.
void WriteH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder, const std::string outputName){
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
	if(outputName.empty()){
		filename = "virasoro_" + to_string(c, 3) + "_" + to_string(hl, 3) + "_" + to_string(hh, 3) + "_" + to_string(hp, 1) + "_" + std::to_string(maxOrder) + ".txt";
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
