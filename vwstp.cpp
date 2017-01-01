#include <string>
#include <vector>
#include "wstp.h"
#include "runfile.h"

extern "C"{
	extern void RunFromComponents(const char* c, const char* hl, const char* hh, const char* hp, const char* maxOrder);
	extern void RunFromFile(const char* filename);
}

void Startup(Runfile_c runfile);

mpc_rnd_t mpfc_class::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpfc_class::default_prec = 64;
mpf_class mpfc_class::tolerance = 1e-10;

int maxThreads = 8;
int precision = 768;
mpf_class tolerance = 1e-10;
bool showProgressBar = false;

int main(int argc, char* argv[]){
	mpf_set_default_prec(precision);
	mpfc_class::set_default_prec(precision);
	mpfc_class::set_default_rnd_mode(MPC_RNDZZ);
	mpfc_class::set_tolerance(tolerance);
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
