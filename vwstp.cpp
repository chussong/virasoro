#include <string>
#include <vector>
#include <fstream>
#include "wstp.h"
#include "runfile.h"

extern "C"{
	extern void RunFromComponents(const char* c, const char* hl, const char* hh, const char* hp, const char* maxOrder);
	extern void RunFromFile(const char* filename);
}

void Startup(Runfile_c runfile);
void ReadDefaults(std::string filename);

mpc_rnd_t mpcomplex::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpcomplex::default_prec = 64;
mpfr::mpreal mpcomplex::tolerance = 1e-10;

int maxThreads = 8;
int precision = 768;
mpfr::mpreal tolerance = 1e-10;
bool showProgressBar = false;

int main(int argc, char* argv[]){
	ReadDefaults("config.txt");
	mpfr::mpreal::set_default_prec(precision);
	mpfr::mpreal::set_default_rnd(MPFR_RNDZ);
	mpcomplex::set_default_prec(precision);
	mpcomplex::set_default_rnd_mode(MPC_RNDZZ);
	mpcomplex::set_tolerance(tolerance);
	return WSMain(argc, argv);
}

void RunFromComponents(const char* c, const char* hl, const char* hh, const char* hp, const char* maxOrder){
	std::vector<std::string> argVec;
	argVec.emplace_back(c);
	argVec.emplace_back(hl);
	argVec.emplace_back(hh);
	argVec.emplace_back(hp);
	argVec.emplace_back(maxOrder);
	Runfile_c runfile(argVec);
	Startup(runfile);
	return;
}

void RunFromFile(const char* filename){
	Runfile_c runfile(filename);
	Startup(runfile);
	return;
}

void ReadDefaults(std::string filename){
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	if((inStream.rdstate() & std::ifstream::failbit) != 0){
		return;
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
			continue;
		}
	}
	inStream.close();
	return;
}

void Startup(Runfile_c runfile){
	runfile.ReadRunfile();
	runfile.SetMaxThreads(maxThreads);
	runfile.SetPrecision(precision);
	runfile.SetTolerance(tolerance);
	runfile.SetProgressBar(showProgressBar);
	WSPutFunction(stdlink, "Map", 2);
	WSPutSymbol(stdlink, "ToExpression");
	WSPutFunction(stdlink, "List", 2*runfile.NumberOfRuns());
	runfile.Execute("w");
	return;
}
