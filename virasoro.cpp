#include "virasoro.h"

// namespace virasoro {

int maxThreads;
int precision;
mpfr::mpreal tolerance;
bool showProgressBar;

constexpr int DEFAULT_THREADS = 8;
constexpr int DEFAULT_PREC = 768;
constexpr double DEFAULT_TOLERANCE = 1e-100;

int core(int argc, char** argv, const bool wolfram){
	ReadDefaults(std::string(getenv("HOME"))+"/.config/virasoro_defaults.txt", wolfram);
	std::vector<std::string> args = CollectArgs(argc, argv);
	std::string options = ParseOptions(args);
	if(wolfram) options.append("w");
	mpfr::mpreal::set_default_prec(precision);
	mpfr::mpreal::set_default_rnd(MPFR_RNDZ);
	int exitCode;
	Runfile_c runfile;
	if(args.size() == 1){
		runfile = args[0];
	} else {
		runfile = args;
	}
	runfile.SetMaxThreads(maxThreads);
	runfile.SetPrecision(precision);
	runfile.SetTolerance(tolerance);
	runfile.SetProgressBar(showProgressBar);
	runfile.ReadRunfile();
#ifdef HAVE_WSTP_
	if(wolfram){
		WSPutFunction(stdlink, "Map", 2);
		WSPutSymbol(stdlink, "ToExpression");
		WSPutFunction(stdlink, "List", 2*runfile.NumberOfRuns());
	}
#endif
	exitCode = runfile.Execute(options);
	return exitCode;
}

std::vector<std::string> CollectArgs(int argc, char** argv){
	std::vector<std::string> args;
	for(int i = 1; i <= argc-1; ++i){
		args.emplace_back(argv[i]);
	}
	return args;
}

void ReadDefaults(const std::string filename, const bool quiet){
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	if((inStream.rdstate() & std::ifstream::failbit) != 0){
		if(!quiet)std::cout << "Creating configuration file at ~/.config/virasoro_defaults.txt" << std::endl;
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
	if(!quiet)std::cout << "Reading default configuration from ~/.config/virasoro_defaults.txt" << std::endl;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		if(lines[i-1].size() >= 12 && lines[i-1].substr(0,10).compare("maxThreads") == 0){
			maxThreads = std::stoi(lines[i-1].substr(11));
			if(maxThreads <= 0) maxThreads = DEFAULT_THREADS;
			continue;
		}
		if(lines[i-1].size() >= 11 && lines[i-1].substr(0,9).compare("precision") == 0){
			precision = std::stoi(lines[i-1].substr(10));
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
	inStream.close();
	return;
}

void CreateConfigFile(const std::string filename){
	std::ofstream outStream;
	outStream.open(filename);
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

std::string to_string(const mpfr::mpreal N, int digits){
	if(digits <= 0) return N.toString(-1, 10, MPFR_RNDN);
	std::string output = N.toString(digits, 10, MPFR_RNDN);
	std::size_t ePos = output.find('e');
	if(ePos < std::string::npos) output = output.replace(ePos, 1, "*10^");
	return output;
}

std::string to_string(const std::complex<mpfr::mpreal> N, int digits, int base){
	char* cstr = mpc_get_str(base, std::max(digits,0), N.mpc_srcptr(), MPC_RNDNN);
	std::string output(cstr);
	mpc_free_str(cstr);
	if(digits < 0){
		return output;
	}
	size_t splitLoc = output.find(" ");
	std::string halves[2];
	halves[0] = output.substr(1, splitLoc-1);
	halves[1] = output.substr(splitLoc+1, output.size()-splitLoc-2);
	mpfr::mpreal mpfHalf;
	for(int i = 1; i <= 2; ++i){
		if(halves[i-1] == "+0" || halves[i-1] == "-0"){
			halves[i-1].clear();
		} else {
			if(halves[i-1][0] == '-'){
				mpfHalf = halves[i-1].substr(1);
			} else {
				mpfHalf = halves[i-1];
			}
			if(mpfHalf < STATICTOLERANCE) halves[i-1].clear();
		}
	}
	size_t eLoc, eEnd;
	for(int i = 1; i <= 2; ++i){
		if(halves[i-1].empty()) continue;
		if((eLoc=halves[i-1].find("e")) < std::string::npos){
			eEnd = halves[i-1].find_first_of(" )", eLoc+3);
			int exp = std::stoi(halves[i-1].substr(eLoc+1, eEnd-eLoc-1));
			if(digits == 0 || std::abs(exp) < digits){
				halves[i-1].erase(eLoc, eEnd-eLoc);
				if(exp > 0){
					if(halves[i-1][0] == '-'){
						halves[i-1].erase(2, 1);
						halves[i-1].insert(exp+2, ".");
					} else {
						halves[i-1].erase(1, 1);
						halves[i-1].insert(exp+1, ".");
					}
				} else {
					if(halves[i-1][0] == '-'){
						halves[i-1].erase(2, 1);
						halves[i-1].insert(1, "0.");
						halves[i-1].insert(3, -1-exp, '0');
					} else {
						halves[i-1].erase(1, 1);
						halves[i-1].insert(0, "0.");
						halves[i-1].insert(2, -1-exp, '0');
					}
				}
			} else {
				if(halves[i-1][eLoc+1] == '+'){
					halves[i-1].replace(eLoc, 2, "*10^");
				} else {
					halves[i-1].replace(eLoc, 1, "*10^");
				}
			}
		}
		eEnd = halves[i-1].find_last_not_of("0");
		if(halves[i-1][eEnd] == '.') halves[i-1].erase(eEnd);
	}
	if(halves[0].empty()){
		if(halves[1].empty()){
			return "0";
		} else {
			halves[1].append("*I");
			return halves[1];
		}
	} else {
		if(halves[1].empty()){
			return halves[0];
		} else {
			if(halves[1][0] == '-') return halves[0] + "-" + halves[1].substr(1) + "*I";
			return halves[0] + "+" + halves[1] + "*I";
		}
	}
}
// } // namespace virasoro
