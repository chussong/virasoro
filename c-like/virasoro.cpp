#include "virasoro.h"

typedef std::chrono::high_resolution_clock Clock;

const int maxThreads = 8;
const int precision = 512;
mpf_t* powOverflow;

int main(int argc, char** argv){
	mpf_set_default_prec(precision);
	auto programStart = Clock::now();
	switch(argc){
		default:	printf("Error: input either the name of a runfile or the five parameters c, hl, hh, hp, maxOrder\n");
					return EXIT_FAILURE;
		case 2:		{mpf_t** runs = (mpf_t**)malloc(sizeof(*runs));
					int* maxOrders = (int*)malloc(sizeof(*maxOrders));
					int numberOfRuns = ReadRunfile(argv[1], runs, maxOrders);
					if(numberOfRuns == -1){
						perror("Error: if you give one argument it must be the name of a runfile.\n");
						return EXIT_FAILURE;
					}
					if(numberOfRuns == 0){
						perror("Error: format of given runfile is invalid.\n");
						return EXIT_FAILURE;
					}
					for(int i = 1; i <= numberOfRuns; ++i){
						std::cout << "Run " << i << ": ";
						for(int j = 1; j <= 4; ++j){
							std::cout << to_string(runs[i-1][j-1], 0) << " ";
						}
						std::cout << maxOrders[i-1];
						std::cout << std::endl;
					}
					int highestMax = 0;
					for(int i = 1; i <= numberOfRuns; ++i){
						if(maxOrders[i-1] > highestMax) highestMax = maxOrders[i-1];
					}
					SetPowOverflow(highestMax);
					auto runStart = Clock::now();
					for(int run = 1; run <= numberOfRuns; ++run){
						runStart = Clock::now();
//						DebugPrintRunVector(runs[run-1], maxOrders[run-1]);
						printf("Beginning run %i.\n", run);
						FindCoefficients(runs[run-1], maxOrders[run-1], true);
						ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);						
					}
					}break;
		case 6:		{mpf_t runVector[4];
					for(int i = 1; i <= 4; ++i){
						mpf_init_set_str(runVector[i-1], argv[i], 10);
					}
					unsigned short int maxOrder = std::atoi(argv[5]);
					maxOrder -= (maxOrder % 2);
					SetPowOverflow(maxOrder);
//					DebugPrintRunVector(runVector, maxOrder);
					FindCoefficients(runVector, maxOrder, false);
					}break;
	}

	ShowTime("Entire computation", programStart);
	return EXIT_SUCCESS;
}

int ReadRunfile(char* filename, mpf_t** &runs, int* &maxOrders){
	int numberOfLines = 0;
	int readStatus = 0;
	FILE* runfile = fopen(filename, "r");
	if(runfile == NULL){

		return -1;
	}
	int c = fgetc(runfile);
	while(c != EOF){
		if(c == '\n') ++numberOfLines;
		c = fgetc(runfile);
	}
	rewind(runfile);
	runs = new mpf_t*[numberOfLines];
	maxOrders = new int[numberOfLines];
	for(int currentLine = 1; currentLine <= numberOfLines; ++currentLine){
		runs[currentLine-1] = new mpf_t[5];
		for(int i = 1; i <= 4; ++i){
			readStatus = ReadMPF(runs[currentLine-1][i-1], runfile);
			if(readStatus != 0) return 0;
		}
		if(readStatus != 0) return 0;
		maxOrders[currentLine-1] = ReadMaxOrder(runfile);
		if(maxOrders[currentLine-1] <= 0) return 0;
	}
	return numberOfLines;
}

int ReadMPF(mpf_t& output, FILE* runfile){
	char inputBuffer[256]{0};
	int c = fgetc(runfile);
	int place = 0;
	while((c >= '0' && c <= '9') || c == '.'){
		inputBuffer[place] = (char)c;
		c = fgetc(runfile);
		++place;
	}
	if(place == 0) return -1;
	if(mpf_init_set_str(output, inputBuffer, 10) == -1) return -1;
	switch(c){
		case ' ':	return 0;
					break;
		case ',':	return 0;
					break;
		case ';':	return 0;
					break;
		case '\n':	return 1;
					break;	
		case EOF:	return 2;
					break;
	}
	return -1;
}

int ReadMaxOrder(FILE* runfile){
	int maxOrder = 0;
	int c = fgetc(runfile);
	while(c >= '0' && c <= '9'){
		maxOrder *= 10;
		maxOrder += c - '0';
		c = fgetc(runfile);
	}
	maxOrder -= maxOrder%2;
	return maxOrder;
}

void SetPowOverflow(unsigned short int maxOrder){
	powOverflow = new mpf_t[maxOrder/256+1];
	mpf_init_set_ui(powOverflow[0], 1);
	for(int i = 1; i <= maxOrder/256; ++i){
		mpf_init_set_ui(powOverflow[i], 16);
		mpf_pow_ui(powOverflow[i], powOverflow[i], 256*i);
	}
}

void DebugPrintRunVector(const mpf_t* runVector, const unsigned short int maxOrder){
	std::cout << "If the code were running right now, it would be initiating a run with c = " << to_string(runVector[0], 0) << ", hl = " << to_string(runVector[1], 0) << ", hh = " << to_string(runVector[2], 0) << ", hp = " << to_string(runVector[3], 0) << ", maxOrder = " << maxOrder << "." << std::endl;
}

void FindCoefficients(const mpf_t* runVector, const unsigned short int maxOrder, const bool multirun){
//	auto timeStart = Clock::now();
	int mnLocation[maxOrder]; 	/* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	int mnMultiplicity[maxOrder] = {0};	/* number of mn combinations giving i+1 = m*n */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int mTable[numberOfMN] = {0};
	unsigned short int nTable[numberOfMN] = {0};
	int mnLookup[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, numberOfMN, mnLocation, mnMultiplicity, maxOrder);
	
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	mpf_t bsq, invBsq, llsq, lhsq, temp1, temp2;
	mpf_inits(bsq, invBsq, llsq, lhsq, temp1, temp2, NULL);
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);
	
	Cpqmn_t Cpqmn(&bsq, &invBsq, numberOfMN, maxOrder, mTable, nTable, mnLookup);
	Cpqmn.FillHpmn();
	
	auto time1 = Clock::now();
	Cpqmn.FillAmn();
	auto time2 = Clock::now();
//	std::cout << "Amn filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// combine Amn into Rmn
	time1 = Clock::now();
	Cpqmn.FillRmn(&llsq, &lhsq);
	time2 = Clock::now();
//	std::cout << "Rmn filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// combine Rmn and hpmn into computation of H
	Hmn_t Hmn(&Cpqmn, numberOfMN, maxOrder, mnLocation, mnMultiplicity, mnLookup);
	time1 = Clock::now();
	Hmn.FillHmn();
	time2 = Clock::now();
//	std::cout << "Hmn filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// corral H by q order and display coefficients
	mpf_t H[maxOrder/2+1];
	mpf_init_set_ui(H[0], 1);
	for(int i = 1; i <= maxOrder/2; ++i) mpf_init(H[i]);
	FillH(H, &Hmn, &Cpqmn, runVector[3], mnLocation, mnMultiplicity, maxOrder);
	if(!multirun) DisplayH(H, runVector[0], runVector[1], runVector[2], runVector[3], maxOrder);
	WriteH(H, runVector[0], runVector[1], runVector[2], runVector[3], maxOrder, multirun);
//	ShowTime(std::string("Computing maxOrder = ").append(std::to_string(maxOrder)), timeStart);
	mpf_clears(bsq, invBsq, llsq, lhsq, temp1, temp2, NULL);
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

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
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
		for(int n = 1; (m+1)*n <= maxOrder; ++n){	// even m
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

void ConvertInputs(mpf_t& bsq, mpf_t& invBsq, mpf_t& llsq, mpf_t& lhsq, const mpf_t& c, const mpf_t& hl, const mpf_t& hh, mpf_t& temp1, mpf_t& temp2){
	mpf_mul(temp1, c, c);
	mpf_mul_ui(temp2, c, 26);
	mpf_sub(temp1, temp1, temp2);
	mpf_add_ui(temp1, temp1, 25);
	mpf_sqrt(temp1, temp1);
	mpf_sub(temp1, c, temp1);
	mpf_sub_ui(temp1, temp1, 13);
	mpf_div_ui(bsq, temp1, 12);
	mpf_ui_div(invBsq, 1, bsq);

	mpf_sub_ui(temp1, c, 1);
	mpf_div_ui(temp1, temp1, 24);
	mpf_sub(llsq, hl, temp1);
	mpf_sub_ui(temp1, c, 1);
	mpf_div_ui(temp1, temp1, 24);
	mpf_sub(lhsq, hh, temp1);
}

void FillH(mpf_t* H, const Hmn_t* Hmn, const Cpqmn_t* Cpqmn, const mpf_t hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	mpf_t temp1, temp2;
	mpf_inits(temp1, temp2, NULL);
	for(int order = 2; order <= maxOrder; order+=2){
		for(int power = 2; power <= order; power+=2){
			for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
				mpf_sub(temp1, hp, Cpqmn->hpmn[scanPos-1]);
				mpf_div(temp1, Cpqmn->Rmn[scanPos-1], temp1);
				mpf_mul(temp1, temp1, Hmn->Hmn[(order-power)/2][scanPos-1]);
				mpf_mul(temp1, temp1, powOverflow[power/256]);
				mpf_set_ui(temp2, 16);
				mpf_pow_ui(temp2, temp2, power%256);
				mpf_mul(temp1, temp1, temp2);
				mpf_add(H[order/2], H[order/2], temp1);
			}
		}
	}
	mpf_clears(temp1, temp2, NULL);
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

std::string to_string(const mpf_t N, int digits){
	if(digits < 0) digits = -digits;	
	char* buffer = new char;
	mp_exp_t* dotPos = new long;
	buffer = mpf_get_str(NULL, dotPos, 10, digits, N);
	std::string output = std::string(buffer);
	if(output.empty()) output.append("0");
	while((long)*dotPos > (int)output.size()) output.append("0");
	if(digits == 0 && (long)*dotPos > 0 && (long)*dotPos < (int)output.size()){
		if(mpf_sgn(N) == -1) *dotPos += 1;
		output.insert((long)*dotPos, ".");
	}
	if(digits > 0 && (long)*dotPos > digits){
		while((int)output.size() < digits) output.append("0");
		if(mpf_sgn(N) == 1) output.insert(1, ".");		
		if(mpf_sgn(N) == -1) output.insert(2, ".");
		output.append("*10^");
		sprintf(buffer, "%ld", (long)*dotPos);
		output.append(buffer);
	}
	delete buffer;
	delete dotPos;
	return output;
}

void DisplayH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder){
	std::cout << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
	}
}

void WriteH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder, const bool multirun){
	std::ofstream outputFile;
	std::string filename = "virasoro_" + to_string(c, 3) + "_" + to_string(hl, 3) + "_" + to_string(hh, 3) + "_" + to_string(hp, 1) + "_" + std::to_string(maxOrder) + ".txt";
	outputFile.open (filename);
	if(!multirun) outputFile << "Given the parameters" << std::endl;
	outputFile << "c = " << to_string(c, 0) << ", h_L = " << to_string(hl, 0) << ", h_H = " << to_string(hh, 0) << ", h_p = " << to_string(hp, 0) << std::endl;	
	if(!multirun){
	outputFile << "the Virasoro block coefficients are as follows:" << std::endl;
		for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
			outputFile << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
		}
	}
	outputFile << "{1";
	for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
		outputFile << "," << to_string(H[orderBy2], 0);
	}
	outputFile << "}" << std::endl;
	outputFile.close();
}
