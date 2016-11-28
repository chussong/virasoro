#include "virasoro.h"

int maxThreads = 8;			// Maximum number of simultaneous threads
int precision = 512;				// Precision of mpf_class in bits
const mpf_class tolerance(1e-20);	// Smaller than this is taken to be 0 for comparisons

mpc_rnd_t mpfc_class::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpfc_class::default_prec = 64;

int main(int argc, char** argv){
	auto programStart = Clock::now();
	std::string options = ParseOptions(argc, argv);
	mpf_set_default_prec(precision);
	mpfc_class::set_default_prec(precision);
	mpfc_class::set_default_rnd_mode(MPC_RNDZZ);
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
		} else if(strcmp(argv[i], "-bb") == 0){
			options.append("bb");
			--newArgc;
			realArg[i] = false;
		} else if(strcmp(argv[i], "-b") == 0){
			options.append("b");
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
				std::cout << to_string(runfile.runs[i-1][j-1], 4) << " ";
			}
			if(runfile.runs[i-1].size() > 4) std::cout << "{";
			for(unsigned int j = 4; j <= runfile.runs[i-1].size(); ++j){
				std::cout << to_string(runfile.runs[i-1][j-1], 4) << ",";
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
	auto runStart = Clock::now();
	std::string outputName;
	std::vector<mpfc_class> complexRunVector;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
		for(unsigned int run = 1; run < runfile.runs.size(); ++run){
			if(bGiven == 0 && runfile.runs[run-1][0] < 25 && runfile.runs[run-1][0] > 1){
				for(int i = 1; i <= 4; ++i) complexRunVector.emplace_back(runfile.runs[run-1][i-1]);
				FindCoefficients<mpfc_class>(complexRunVector, runfile.maxOrders[run-1], outputName, bGiven);
			} else {
				FindCoefficients<mpf_class>(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			}
			std::cout << ",";
		}
		if(bGiven == 0 && runfile.runs[runfile.runs.size()-1][0] > 1 && runfile.runs[runfile.runs.size()-1][0] > 1){
			for(int i = 1; i <= 4; ++i) complexRunVector.emplace_back(runfile.runs[runfile.runs.size()-1][i-1]);
			FindCoefficients<mpfc_class>(complexRunVector, runfile.maxOrders[runfile.runs.size()-1], outputName, bGiven);
		} else {
			FindCoefficients<mpf_class>(runfile.runs[runfile.runs.size()-1], runfile.maxOrders[runfile.runs.size()-1], outputName, bGiven);
		}
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
			if(bGiven == 0 && runfile.runs[run-1][0] < 25 && runfile.runs[run-1][0] > 1){
				for(int i = 1; i <= 4; ++i) complexRunVector.emplace_back(runfile.runs[run-1][i-1]);
				FindCoefficients<mpfc_class>(complexRunVector, runfile.maxOrders[run-1], outputName, bGiven);
			} else {
				FindCoefficients<mpf_class>(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			}
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
	if(bGiven == 0 && runVector[0] < 25 && runVector[0] > 1){
		std::vector<mpfc_class> complexRunVector;
		for(int i = 1; i <= 4; ++i) complexRunVector.emplace_back(runVector[i-1]);
		FindCoefficients<mpfc_class>(complexRunVector, maxOrder, outputName, bGiven);
	} else {
		FindCoefficients<mpf_class>(runVector, maxOrder, outputName, bGiven);
	}
	if(wolframOutput) std::cout << "}";
	return 0;
}

void DebugPrintRunVector(const mpf_class* runVector, const std::vector<mpf_class> hp, const unsigned short int maxOrder){
	for(unsigned int i = 1; i <= hp.size(); ++i){
		std::cout << "Pretend this is a run with c = " << to_string(runVector[0], 0) << ", hl = " << to_string(runVector[1], 0) << ", hh = " << to_string(runVector[2], 0) << ", hp = " << to_string(hp[i-1], 0) << ", maxOrder = " << maxOrder << "." << std::endl;
	}
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
	double Nd = N.get_d();
	if(output.empty()) return "0";
	if(digits > 0 && dotPos < digits){				// number small enough, just write it
		while(dotPos > (int)output.size()) output.append("0");
		if(dotPos < (int)output.size()){
			if(Nd >= 1){
				output.insert(abs(dotPos), ".");
			} else if(Nd < 1 && Nd > 0) {
				output.insert(0, "0.");
			} else if(Nd < 0) {
				output.insert(abs(dotPos)+1, ".");
			}
		}
	}
	if(digits > 0 && abs(dotPos) > digits){			// number too big, use a*10^b
		while((int)output.size() < digits) output.append("0");
		if(sgn(N) == 1) output.insert(1, ".");		
		if(sgn(N) == -1) output.insert(2, ".");
		output.append("*10^");
		output.append(std::to_string(dotPos));
	}
	if(digits == 0){								// entire number has been requested
		if(Nd >= 1){
			while(dotPos > (int)output.size()) output.append("0");
			output.insert(dotPos, ".");
		} else if(Nd > 0) {
			while(-dotPos > (int)output.size()) output.insert(0, "0");
			output.insert(0, "0.");
		} else if(Nd > -1) {
			while(-dotPos > (int)output.size()) output.insert(1, "0");
			output.insert(1, "0.");
		} else if(Nd <= -1) {
			while(dotPos > (int)output.size()) output.append("0");
			output.insert(dotPos+1, ".");
		}
	}
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
