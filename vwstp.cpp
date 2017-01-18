#include <cstdio>
#include "wstp.h"
#include "virasoro.h"

extern "C"{
	extern void RunFromComponents(const char* c, const char* hl, const char* hh, const char* hp, const char* maxOrder);
	extern void RunFromFile(const char* filename);
}

int main(int argc, char* argv[]){
	return WSMain(argc, argv);
}

void RunFromComponents(const char* c, const char* hl, const char* hh, const char* hp, const char* maxOrder){
/*	std::vector<std::string> argVec;
	argVec.emplace_back(c);
	argVec.emplace_back(hl);
	argVec.emplace_back(hh);
	argVec.emplace_back(hp);
	argVec.emplace_back(maxOrder);
	Runfile_c runfile(argVec);
	Startup(runfile);*/
	char* argv[6];
	argv[0] = new char[17];
	strcpy(argv[0], "vwstp_components");
	strcpy(argv[1], c);
	strcpy(argv[2], hl);
	strcpy(argv[3], hh);
	strcpy(argv[4], hp);
	strcpy(argv[5], maxOrder);
	core(6, argv, true);
	delete[] argv[0];
	return;
}

void RunFromFile(const char* filename){
	char* argv[2];
	argv[0] = new char[11];
	strcpy(argv[0], "vwstp_file");
	strcpy(argv[1], filename);
	core(2, argv, true);
	delete[] argv[0];
	return;
}
