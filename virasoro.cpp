#include "virasoro.h"

int maxThreads = 8;					// Maximum number of simultaneous threads
int precision = 512;				// Precision of mpf_class and mpfc_class in bits
mpf_class tolerance(1e-20);			// Smaller than this is taken to be 0 for comparisons
bool showProgressBar = true;		// Show progress bar during FillHmn()

mpc_rnd_t mpfc_class::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpfc_class::default_prec = 64;

int main(int argc, char** argv){
	auto programStart = Clock::now();
	ReadDefaults("config.txt");
	std::vector<std::string> args = CollectArgs(argc, argv);
	std::string options = ParseOptions(args);
	mpf_set_default_prec(precision);
	mpfc_class::set_default_prec(precision);
	mpfc_class::set_default_rnd_mode(MPC_RNDZZ);
	int exitCode;
	switch(args.size()){
		default:	//printf("Error: input either the name of a runfile or the five parameters c, hl, hh, hp, maxOrder\n");
					{Runfile_c runfile(args);
					exitCode = ExecuteRunfile(runfile, options);
					break;}
					//return EXIT_FAILURE;
		case 1:		{Runfile_c runfile(args[0]);
					exitCode = ExecuteRunfile(runfile, options);
					break;}
		case 5:		{Runfile_c runfile(args);
					exitCode = ExecuteRunfile(runfile, options);
					break;}
	}
	if(options.find("m", 0) != std::string::npos) ShowTime("Entire computation", programStart);
	return exitCode;
}

std::vector<std::string> CollectArgs(int argc, char** argv){
	std::vector<std::string> args;
	for(int i = 1; i <= argc-1; ++i){
		args.emplace_back(argv[i]);
	}
	return args;
}

void ReadDefaults(std::string filename){
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	if((inStream.rdstate() & std::ifstream::failbit) != 0){
		CreateConfigFile(filename);
		inStream.open(filename, std::ifstream::in);
	}
	std::string currentLine;
	std::vector<std::string> lines;
	while(true){
		std::getline(inStream, currentLine);
		if(currentLine.empty()) break;
		lines.push_back(currentLine);
	}
	for(unsigned int i = 1; i <= lines.size(); ++i){
		if(lines[i-1].size() >= 12 && lines[i-1].substr(0,10).compare("maxThreads") == 0){
			maxThreads = std::stoi(lines[i-1].substr(11));
			continue;
		}
		if(lines[i-1].size() >= 11 && lines[i-1].substr(0,9).compare("precision") == 0){
			precision = std::stoi(lines[i-1].substr(10));
			continue;
		}
		if(lines[i-1].size() >= 11 && lines[i-1].substr(0,9).compare("tolerance") == 0){
			tolerance = lines[i-1].substr(10).c_str();
			continue;
		}
		if(lines[i-1].size() >= 17 && lines[i-1].substr(0,15).compare("showProgressBar") == 0){
			if(lines[i-1].substr(16) == "false"){
				showProgressBar = false;
			} else {
				showProgressBar = true;
			}
		}
	}
	inStream.close();
	return;
}

void CreateConfigFile(std::string filename){
	std::ofstream outStream;
	outStream.open(filename, std::ofstream::out);
	outStream << "[default parameters]" << std::endl;
	outStream << "maxThreads=8" << std::endl;
	outStream << "precision=512" << std::endl;
	outStream << "tolerance=1e-20" << std::endl;
	outStream << "showProgressBar=true" << std::endl;
	outStream.close();
	return;
}

std::string ParseOptions(std::vector<std::string> &args){
	std::string options = "";
	std::vector<bool> realArg;
	realArg.resize(args.size());
	for(unsigned int i = 1; i <= args.size(); ++i){
		if(args[i-1][0] != '-'){
			realArg[i-1] = true;
			continue;
		} else if(args[i-1].substr(0,2).compare("-m") == 0){
			options.append("m");
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,2).compare("-c") == 0){
			options.append("c");
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,3).compare("-bb") == 0){
			options.append("bb");
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,2).compare("-b") == 0){
			options.append("b");
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,2).compare("-p") == 0){
			precision = std::stoi(args[i-1].substr(2));
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,2).compare("-t") == 0){
			maxThreads = std::stoi(args[i-1].substr(2));
			realArg[i-1] = false;
		} else {
			realArg[i-1] = true;
		}
	}
	for(unsigned int i = realArg.size(); i >= 1; --i){
		if(!realArg[i-1]){
			args.erase(args.begin()+i-1);
		}
	}
	return options;
}

int RunFromFile(std::string filename, const std::string options){
	Runfile_c runfile(filename);
	return ExecuteRunfile(runfile, options);
}

int RunFromTerminal(std::vector<std::string> args, const std::string options){
	Runfile_c runfile(args);
	return ExecuteRunfile(runfile, options);
/*	std::vector<mpf_class> runVector;
	for(int i = 1; i <= 4; ++i){
		runVector.emplace_back(args[i-1]);
	}
	unsigned short int maxOrder = std::stoi(args[4]);
	maxOrder -= (maxOrder % 2);
	DebugPrintRunVector(runVector, hp, maxOrder);
	std::string outputName;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
	} else if(consoleOutput) {
		outputName = "__CONSOLE";
	} else {
		outputName = NameOutputFile("");
		std::remove(outputName.c_str());
	}
	if(bGiven == 0 && runVector[0] < 25 && runVector[0] > 1){
		std::vector<mpfc_class> complexRunVector;
		for(unsigned int i = 1; i <= runVector.size(); ++i) complexRunVector.emplace_back(runVector[i-1]);
		FindCoefficients<mpfc_class>(complexRunVector, maxOrder, outputName, bGiven);
		complexRunVector.clear();
	} else {
		FindCoefficients<mpf_class>(runVector, maxOrder, outputName, bGiven);
	}
	if(wolframOutput) std::cout << "}";
	return 0;*/
}

int ExecuteRunfile(Runfile_c runfile, std::string options){
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;	
	int bGiven = 0;
	if(options.find("b", 0) != std::string::npos) bGiven = 1;
	if(options.find("bb", 0) != std::string::npos) bGiven = 2;
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
		std::cout << "Output will be saved to " << runfile.NameOutputFile() << ". If it exists, it will be overwritten." << std::endl;
	}
	int highestMax = 0;
	for(unsigned int i = 1; i <= runfile.runs.size(); ++i){
		if(runfile.maxOrders[i-1] > highestMax) highestMax = runfile.maxOrders[i-1];
	}
	auto runStart = Clock::now();
	std::string outputName;
	std::vector<mpf_class> realRunVector;
	bool allReal;
	if(wolframOutput){
		showProgressBar = false;
		outputName = "__MATHEMATICA";
		std::cout << "{";
		for(unsigned int run = 1; run <= runfile.runs.size(); ++run){
			allReal = true;
			for(unsigned int i = 1; i <= runfile.runs[run-1].size(); ++i){
				if(!runfile.runs[run-1][i-1].isReal()){
					allReal = false;
					break;
				}
			}
			if(bGiven == 0 && runfile.runs[run-1][0].realPart() < 25 && runfile.runs[run-1][0].realPart() > 1) allReal = false;
			if(allReal){
				for(unsigned int i = 1; i <= runfile.runs[run-1].size(); ++i) realRunVector.push_back(runfile.runs[run-1][i-1].realPart());
				FindCoefficients<mpf_class>(realRunVector, runfile.maxOrders[run-1], outputName, bGiven);
				realRunVector.clear();
			} else {
				FindCoefficients<mpfc_class>(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			}
			if(run < runfile.runs.size()) std::cout << ",";
		}
		std::cout << "}";
	} else {
		if(consoleOutput){
			outputName = "__CONSOLE";
		} else {
			outputName = runfile.NameOutputFile();
			std::remove(outputName.c_str());
		}
		for(unsigned int run = 1; run <= runfile.runs.size(); ++run){
			runStart = Clock::now();
//			DebugPrintRunVector(runs[run-1], hp[run-1], maxOrders[run-1]);
			std::cout << "Beginning run " << run << " of " << runfile.runs.size() << "." << std::endl;
			allReal = true;
			for(unsigned int i = 1; i <= runfile.runs[run-1].size(); ++i){
				if(!runfile.runs[run-1][i-1].isReal()){
					allReal = false;
					break;
				}
			}
			if(bGiven == 0 && runfile.runs[run-1][0].realPart() < 25 && runfile.runs[run-1][0].realPart() > 1) allReal = false;
			if(allReal){
				for(unsigned int i = 1; i <= runfile.runs[run-1].size(); ++i) realRunVector.push_back(runfile.runs[run-1][i-1].realPart());
				FindCoefficients<mpf_class>(realRunVector, runfile.maxOrders[run-1], outputName, bGiven);
				realRunVector.clear();
			} else {
				FindCoefficients<mpfc_class>(runfile.runs[run-1], runfile.maxOrders[run-1], outputName, bGiven);
			}
			ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);
		}
	}
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
	mp_exp_t expo, dotPos;
	std::string output = N.get_str(expo, 10, digits);
	dotPos = expo;
	--expo;
	if(dotPos >= 0 && sgn(N) == -1) ++dotPos;
	if(dotPos < 0 && sgn(N) == -1) --dotPos;
	double Nd = N.get_d();
	if(output.empty()) return "0";
	if(digits > 0 && abs(expo) <= digits-1){			// number small enough, just write it
		while(dotPos > (int)output.size()) output.append("0");
		if(dotPos < (int)output.size()){
			if(Nd >= 1 || Nd <= -1){
				output.insert(dotPos, ".");
			} else if(Nd < 1 && Nd > 0) {
				output.insert(0, std::max(-(int)dotPos,0), '0');
				output.insert(0, "0.");
			} else if(Nd < 0 && Nd > -1){
				output.insert(1, std::max(-(int)dotPos,0), '0');
				output.insert(1, "0.");
			}
		}
	}
	if(digits > 0 && abs(expo) > digits-1){			// number too big, use a*10^b
		while((int)output.size() < digits) output.append("0");
		if(sgn(N) == 1) output.insert(1, ".");		
		if(sgn(N) == -1){
			output.append("0");
			output.insert(2, ".");
		}
		output.append("*10^");
		output.append(std::to_string(expo));
	}
	if(digits == 0){								// entire number has been requested
		if(Nd >= 1 || Nd <= -1){
			while(dotPos > (int)output.size()) output.append("0");
			if(dotPos < (int)output.size()) output.insert(dotPos, ".");
		} else if(Nd > 0) {
			output.insert(0, std::max(-(int)dotPos, 0), '0');
			output.insert(0, "0.");
		} else if(Nd > -1) {
			output.insert(1, std::max(-(int)dotPos, 0), '0');
			output.insert(1, "0.");
		}
	}
	return output;
}
