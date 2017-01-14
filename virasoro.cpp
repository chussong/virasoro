#include "virasoro.h"
#ifdef HAVE_WSTP_
WSLINK stdlink = 0;
#endif

mpc_rnd_t mpfr::mpcomplex::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpfr::mpcomplex::default_prec = 64;
mpfr::mpreal mpfr::mpcomplex::tolerance = 1e-10;

int maxThreads;
int precision;
mpfr::mpreal tolerance;
bool showProgressBar;

constexpr int DEFAULT_THREADS = 8;
constexpr int DEFAULT_PREC = 768;
constexpr double DEFAULT_TOLERANCE = 1e-10;

int main(int argc, char** argv){
	ReadDefaults(std::string(getenv("HOME"))+"/.config/virasoro_defaults.txt");
	std::vector<std::string> args = CollectArgs(argc, argv);
	std::string options = ParseOptions(args);
	mpfr::mpreal::set_default_prec(precision);
	mpfr::mpreal::set_default_rnd(MPFR_RNDZ);
	mpfr::mpcomplex::set_default_prec(precision);
	mpfr::mpcomplex::set_default_rnd_mode(MPC_RNDZZ);
	mpfr::mpcomplex::set_tolerance(tolerance);
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

void ReadDefaults(std::string filename){
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	if((inStream.rdstate() & std::ifstream::failbit) != 0){
		std::cout << "Creating configuration file at ~/.config/virasoro_defaults.txt" << std::endl;
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
	std::cout << "Reading default configuration from ~/.config/virasoro.defaults.txt" << std::endl;
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

void CreateConfigFile(std::string filename){
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
