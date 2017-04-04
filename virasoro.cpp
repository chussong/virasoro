#include "virasoro.h"

namespace virasoro {

int maxThreads;
int precision;
mpfr::mpreal tolerance;
bool showProgressBar;

constexpr int DEFAULT_THREADS = 8;
constexpr int DEFAULT_PREC = 768;
constexpr double DEFAULT_TOLERANCE = 1e-100;

int core(int argc, char** argv, const bool wolfram){
	bool quiet = wolfram;
	std::vector<std::string> args = CollectArgs(argc, argv);
	std::string options = ParseOptions(args);
	if(options.find('v') != std::string::npos){
		VersionCheck();
		return EXIT_SUCCESS;
	}
	if(wolfram) options.append("w");
	if(options.find('m') != std::string::npos) quiet = true;
	ReadDefaults(std::string(getenv("HOME"))+"/.config/virasoro_defaults.txt", quiet);
	DoOptions(options, quiet);
	mpfr::mpreal::set_default_prec(precision);
	mpfr::mpreal::set_default_rnd(MPFR_RNDZ);
	int exitCode;
	if(args.size() == 1){
		Runfile_c runfile(args[0]);
		runfile.ReadRunfile();
#ifdef HAVE_WSTP
		if(wolfram){
			WSPutFunction(stdlink, "Map", 2);
			WSPutSymbol(stdlink, "ToExpression");
			WSPutFunction(stdlink, "List", 2*runfile.NumberOfRuns());
		}
#endif
		exitCode = DoRuns(runfile, options);
	} else {
		Runfile_c runfile(args);
		runfile.ReadRunfile();
#ifdef HAVE_WSTP
		if(wolfram){
			WSPutFunction(stdlink, "Map", 2);
			WSPutSymbol(stdlink, "ToExpression");
			WSPutFunction(stdlink, "List", 2*runfile.NumberOfRuns());
		}
#endif
		exitCode = DoRuns(runfile, options);
	}
	return exitCode;
}

std::vector<std::string> CollectArgs(int argc, char** argv){
	std::vector<std::string> args;
	for(int i = 1; i <= argc-1; ++i){
		args.emplace_back(argv[i]);
	}
	return args;
}

void UseCompileTimeDefaults(){
	perror("Error: configuration file at ~/.config/virasoro_defaults.txt was "
			"unreadable. Using compile-time defaults.\n");
	maxThreads = DEFAULT_THREADS;
	precision = DEFAULT_PREC;
	tolerance = DEFAULT_TOLERANCE;
	showProgressBar = true;
}

void ReadDefaults(const std::string filename, const bool quiet){
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	if((inStream.rdstate() & std::ifstream::failbit) != 0){
		if(!quiet)std::cout << "Creating configuration file at ~/.config/virasoro_defaults.txt" << std::endl;
		try{
			CreateConfigFile(filename);
		} catch(std::runtime_error) {
			UseCompileTimeDefaults();
			return;
		}
		inStream.open(filename, std::ifstream::in);
	}
	if(!quiet)std::cout << "Reading default configuration from ~/.config/virasoro_defaults.txt" << std::endl;
	std::string currentLine;
	std::vector<std::string> lines;
	while(true){
		std::getline(inStream, currentLine);
		if(currentLine.empty()) break;
		lines.push_back(currentLine);
	}
	inStream.close();
	if(lines.size() == 0){
		UseCompileTimeDefaults();
		return;
	}
	for(unsigned int i = 1; i <= lines.size(); ++i){
		if(lines[i-1].size() >= 12 && lines[i-1].substr(0,10).compare("maxThreads") == 0){
			try{
				maxThreads = std::stoi(lines[i-1].substr(11));
			} catch(std::invalid_argument){
				perror("Error: maxThreads in ~/.config/virasoro_defaults.txt is"
						"not a number. Please follow the syntax in the README.\n");
				throw(std::invalid_argument("maxThreads"));
			}
			if(maxThreads <= 0) maxThreads = DEFAULT_THREADS;
			continue;
		}
		if(lines[i-1].size() >= 11 && lines[i-1].substr(0,9).compare("precision") == 0){
			try{
				precision = std::stoi(lines[i-1].substr(10));
			} catch(std::invalid_argument){
				perror("Error: precision in ~/.config/virasoro_defaults.txt is"
						"not a number. Please follow the syntax in the README.\n");
				throw(std::invalid_argument("precision"));
			}
			if(precision <= 0) precision = DEFAULT_PREC;
			continue;
		}
		if(lines[i-1].size() >= 11 && lines[i-1].substr(0,9).compare("tolerance") == 0){
			tolerance = lines[i-1].substr(10).c_str();
			if(tolerance <= 0) tolerance = DEFAULT_TOLERANCE;
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
	return;
}

void CreateConfigFile(const std::string filename){
	std::ofstream outStream;
	outStream.open(filename);
	if(outStream.fail()){
		perror("Error: failed to create configuration file at "
				"~/.config/virasoro_defaults.txt. Do you have permission to "
				"write to this directory?");
		throw(std::runtime_error("CreateConfigFile"));
	}
	outStream << "[default parameters]" << std::endl;
	outStream << "maxThreads=" << DEFAULT_THREADS << std::endl;
	outStream << "precision=" << DEFAULT_PREC << std::endl;
	outStream << "tolerance=" << DEFAULT_TOLERANCE << std::endl;
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
		} else if(args[i-1].substr(0,2).compare("-v") == 0){
			options.append("v");
			realArg[i-1] = false;
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
			options.append(args[i-1]);
			realArg[i-1] = false;
		} else if(args[i-1].substr(0,2).compare("-t") == 0){
			options.append(args[i-1]);
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

void DoOptions(const std::string& options, const bool quiet){
	std::size_t optionLoc, optionEnd;
	if((optionLoc = options.find('p')) != std::string::npos){
		optionEnd = options.find_first_not_of("0123456789", optionLoc+1) - 1;
		try{
			precision = std::stoi(options.substr(optionLoc+1,optionEnd-optionLoc));
		} catch(std::invalid_argument){
			perror("Error: precision specified incorrectly. Proper syntax for"
					"this option is -p#, e.g -p256.\n");
			throw(std::invalid_argument("precision"));
		}
		if(!quiet) std::cout << "Precision set to " << precision << " bits." << std::endl;
	}
	if((optionLoc = options.find('t')) != std::string::npos){
		optionEnd = options.find_first_not_of("0123456789", optionLoc+1) - 1;
		try{
			maxThreads = std::stoi(options.substr(optionLoc+1,optionEnd-optionLoc));
		} catch(std::invalid_argument){
			perror("Error: maxThreads specified incorrectly. Proper syntax for"
					"this option is -t#, e.g -t8.\n");
			throw(std::invalid_argument("maxThreads"));
		}
		if(!quiet) std::cout << "Number of threads set to " << maxThreads << "." << std::endl;
	}
}

std::string NameOutputFile(const Runfile_c& runfile){
	std::string outputname = runfile.filename;
	if(runfile.lines.size() == 1){
		outputname = "virasoro_" + to_string(runfile.runs[0][0], 3) + "_" + to_string(runfile.runs[0][1], 3) + "_" + to_string(runfile.runs[0][2], 3) + "_" + to_string(runfile.runs[0][3], 1) + "_" + std::to_string(runfile.maxOrders[0]) + ".txt";
	} else {
		std::size_t delPos = outputname.find(".txt");
		if(delPos != std::string::npos) outputname.erase(delPos, 4);
		outputname.append("_results.txt");
	}
	return outputname;
}

bool ParamsReal(const std::vector<std::complex<mpfr::mpreal>>& runVec){
	for(unsigned int i = 0; i < runVec.size(); ++i){
		if(runVec[i].imag() > tolerance) return false;
	}
	return true;
}

void CheckRealityAndRun(const std::vector<std::complex<mpfr::mpreal>>& runVec, const int maxOrder, const std::string outputName, const int bGiven){
	if(ParamsReal(runVec)){
		if(bGiven != 0 || (runVec[0].real() > 25 || runVec[0].real() < 1)){
			std::vector<mpfr::mpreal> realRunVector;
			for(unsigned int i = 0; i < runVec.size(); ++i){
				realRunVector.push_back(runVec[i].real());
			}
			FindCoefficients<mpfr::mpreal>(realRunVector, maxOrder, outputName, bGiven, false);
		} else {
			FindCoefficients<std::complex<mpfr::mpreal>>(runVec, maxOrder, outputName, bGiven, false);
		}
	} else {
		FindCoefficients<std::complex<mpfr::mpreal>>(runVec, maxOrder, outputName, bGiven, true);
	}
}

int DoRuns(const Runfile_c& runfile, const std::string options){
	auto programStart = Clock::now();
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;	
	const bool wstp = options.find("w", 0) != std::string::npos;
	int bGiven = 0;
	if(options.find("b", 0) != std::string::npos) bGiven = 1;
	if(options.find("bb", 0) != std::string::npos) bGiven = 2;
	if(!wstp && !wolframOutput){
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
		if(!consoleOutput) std::cout << "Output will be saved to " << NameOutputFile(runfile) << ". If it exists, it will be overwritten." << std::endl;
	}
	int highestMax = 0;
	for(unsigned int i = 0; i < runfile.maxOrders.size(); ++i){
		if(runfile.maxOrders[i] > highestMax) highestMax = runfile.maxOrders[i];
	}
	auto runStart = Clock::now();
	std::string outputName;
	if(wstp || wolframOutput){
		showProgressBar = false;
		if(wolframOutput){
			outputName = "__MATHEMATICA";
			std::cout << "{";
		} else {
			outputName = "__WSTP";
		}
		for(unsigned int run = 0; run < runfile.runs.size(); ++run){
			CheckRealityAndRun(runfile.runs[run], runfile.maxOrders[run], outputName, bGiven);
			if(wolframOutput && run < runfile.runs.size()-1) std::cout << ",";
		}
		if(wolframOutput) std::cout << "\b}";
	} else {
		if(consoleOutput){
			outputName = "__CONSOLE";
		} else {
			outputName = NameOutputFile(runfile);
			std::remove(outputName.c_str());
		}
		for(unsigned int run = 0; run < runfile.runs.size(); ++run){
			runStart = Clock::now();
			std::cout << "Beginning run " << run+1 << " of " << runfile.runs.size() << "." << std::endl;
			CheckRealityAndRun(runfile.runs[run], runfile.maxOrders[run], outputName, bGiven);

			if(runfile.runs.size() > 1) ShowTime(std::string("Computing run ").append(std::to_string(run+1)), runStart);
		}
	}
	if(!wstp && !wolframOutput) ShowTime("Entire computation", programStart);
	return 0;
}

void ShowTime(const std::string computationName, const std::chrono::time_point<std::chrono::high_resolution_clock> timeStart){
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

void VersionCheck(){
	std::cout << "This is virasoro v" << versionNumber
		<< ", published on " << versionDate << "." << std::endl;
}

} // namespace virasoro
