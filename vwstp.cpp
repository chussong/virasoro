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
	argv[1] = new char[strlen(c)+1];
	strcpy(argv[1], c);
	argv[2] = new char[strlen(hl)+1];
	strcpy(argv[2], hl);
	argv[3] = new char[strlen(hh)+1];
	strcpy(argv[3], hh);
	argv[4] = new char[strlen(hp)+1];
	strcpy(argv[4], hp);
	argv[5] = new char[strlen(maxOrder)+1];
	strcpy(argv[5], maxOrder);
	/*argv[1] = const_cast<char*>(c);
	argv[2] = const_cast<char*>(hl);
	argv[3] = const_cast<char*>(hh);
	argv[4] = const_cast<char*>(hp);
	argv[5] = const_cast<char*>(maxOrder);*/
	virasoro::core(6, argv, true);
	//for(int i = 0; i < 6; ++i) delete[] argv[i];
	return;
}

void RunFromFile(const char* filename){
	char* argv[2];
	argv[0] = new char[11];
	strcpy(argv[0], "vwstp_file");
	argv[1] = new char[strlen(filename)];
	strcpy(argv[1], filename);
	virasoro::core(2, argv, true);
	delete[] argv[0];
	return;
}
