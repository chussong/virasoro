#include "virasoro.h"

typedef std::chrono::high_resolution_clock Clock;

const int maxThreads = 8;
const int precision = 512;
mpf_class* powOverflow;

int main(int argc, char** argv){
	mpf_set_default_prec(precision);
	auto programStart = Clock::now();
	std::string outputName;
	switch(argc){
		default:	printf("Error: input either the name of a runfile or the five parameters c, hl, hh, hp, maxOrder\n");
					return EXIT_FAILURE;
		case 2:		{mpf_class** rawRuns = (mpf_class**)malloc(sizeof(*rawRuns));
					int* rawMaxOrders = (int*)malloc(sizeof(*rawMaxOrders));
					int rawNumberOfRuns = ReadRunfile(argv[1], rawRuns, rawMaxOrders);
					switch(rawNumberOfRuns){
						case 0:		perror("Error: zero valid runs detected in given runfile.\n");
									return EXIT_FAILURE;
						case -1:	perror("Error: if you give one argument it must be the name of a runfile.\n");
									return EXIT_FAILURE;
						case -2:	perror("Error: a curly brace '{' was found indicating a batch run, but it was not followed by a valid macro.\n");
									return EXIT_FAILURE;
						case -3:	perror("Error: expected a number in the runfile but failed to read one.\n");
									return EXIT_FAILURE;
					}
					// Check runs for duplicates, compress runs differing only by hp
					int numberOfRuns = rawNumberOfRuns;
					bool crunched[rawNumberOfRuns]{0};
					std::vector<mpf_class> rawHp[rawNumberOfRuns];
					for(int i = 1; i <= rawNumberOfRuns; ++i){
						rawHp[i-1].push_back(rawRuns[i-1][3]);
						for(int j = i+1; j <= rawNumberOfRuns; ++j){
							if(crunched[j-1]) continue;
							switch(RunCompare(rawRuns[i-1], rawRuns[j-1])){
								case 0:		crunched[j-1] = false;	// runs different, do nothing
											break;				
								case -1:	if(rawMaxOrders[j-1] > rawMaxOrders[i-1]) rawMaxOrders[i-1] = rawMaxOrders[j-1];
											crunched[j-1] = true;	// runs identical, crunch them together
											--numberOfRuns;
											break;
								case -2:	rawHp[i-1].push_back(rawRuns[j-1][3]);	// runs differ by hp
											if(rawMaxOrders[j-1] > rawMaxOrders[i-1]) rawMaxOrders[i-1] = rawMaxOrders[j-1];
											crunched[j-1] = true;
											--numberOfRuns;
											break;
							}
						}
					}
					int runID = 1;
					mpf_class* runs[numberOfRuns];
					std::vector<mpf_class> hp[numberOfRuns];
					int maxOrders[numberOfRuns];
					for(int i = 1; i <= rawNumberOfRuns; ++i){
						if(!crunched[i-1]){
							runs[runID-1] = new mpf_class[3];
							runs[runID-1][0] = rawRuns[i-1][0];
							runs[runID-1][1] = rawRuns[i-1][1];
							runs[runID-1][2] = rawRuns[i-1][2];
							hp[runID-1] = rawHp[i-1];
							maxOrders[runID-1] = rawMaxOrders[i-1];
							delete[] rawRuns[i-1];
							rawHp[i-1].clear();
							++runID;
						}
					}
					std::free(rawRuns);
					std::free(rawMaxOrders);
					for(int i = 1; i <= numberOfRuns; ++i){
						std::cout << "Run " << i << ": ";
						for(int j = 1; j <= 3; ++j){
							std::cout << to_string(runs[i-1][j-1], 0) << " ";
						}
						if(hp[i-1].size() > 1) std::cout << "{";
						for(unsigned int j = 1; j <= hp[i-1].size(); ++j){
							std::cout << to_string(hp[i-1][j-1], 0) << ",";
						}
						if(hp[i-1].size() > 1){
							std::cout << "\b} ";
						} else {
							std::cout << "\b ";
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
					outputName = NameOutputFile(argv[1]);
					std::remove(outputName.c_str());
					for(int run = 1; run <= numberOfRuns; ++run){
						runStart = Clock::now();
//						DebugPrintRunVector(runs[run-1], hp[run-1], maxOrders[run-1]);
						printf("Beginning run %i of %i.\n", run, numberOfRuns);
						FindCoefficients(runs[run-1], hp[run-1], maxOrders[run-1], outputName);
						ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);						
					}
					}break;
		case 6:		{mpf_class runVector[4];
					for(int i = 1; i <= 4; ++i){
						runVector[i-1] = argv[i];
					}
					std::vector<mpf_class> hp;
					hp.push_back(runVector[3]);
					unsigned short int maxOrder = std::atoi(argv[5]);
					maxOrder -= (maxOrder % 2);
					SetPowOverflow(maxOrder);
//					DebugPrintRunVector(runVector, hp, maxOrder);
					outputName = NameOutputFile(nullptr);
					FindCoefficients(runVector, hp, maxOrder, outputName);
					}break;
	}

	ShowTime("Entire computation", programStart);
	return EXIT_SUCCESS;
}

int ReadRunfile(char* filename, mpf_class** &runs, int* &maxOrders){
	int numberOfLines = 0;
	int readStatus = 0;
	int c;
	std::ifstream inStream;
	std::ofstream tempStream, outStream;
	inStream.open(filename, std::ifstream::in);
	outStream.open("__expfile.txt", std::ofstream::out);
	if((inStream.rdstate() & std::ifstream::failbit) != 0 || (outStream.rdstate() & std::ofstream::failbit) != 0){
		return -1;
	}
	std::string currentLine, firstHalf, secondHalf, insideBraces;
	std::size_t lbPos, rbPos, numStart, numEnd;
	mpf_class lowerBound, upperBound, increment;	
	std::getline(inStream, currentLine);
	bool needRerun = false;
	while(true){
		if((lbPos = currentLine.find("{")) != std::string::npos){
			firstHalf = currentLine.substr(0, lbPos);
			if((rbPos = currentLine.find("}")) == std::string::npos){
				return -2;
			} else {
				secondHalf = currentLine.substr(rbPos+1, currentLine.length()-rbPos-1);
				insideBraces = currentLine.substr(lbPos+1, rbPos-lbPos-1);
				numStart = insideBraces.find_first_of("0123456789");
				numEnd = insideBraces.find_first_of(" ,;", numStart+1);
				lowerBound = insideBraces.substr(numStart, numEnd-numStart);
				numStart = insideBraces.find_first_of("0123456789", numEnd+1);
				numEnd = insideBraces.find_first_of(" ,;", numStart+1);
				upperBound = insideBraces.substr(numStart, numEnd-numStart);
				numStart = insideBraces.find_first_of("0123456789", numEnd+1);
				numEnd = insideBraces.find_first_of(" ,;", numStart+1);
				increment = insideBraces.substr(numStart, numEnd-numStart);
				if(lowerBound > upperBound || increment <= 0) return -2;
				for(mpf_class currentValue = lowerBound; currentValue <= upperBound; currentValue += increment){
					outStream << firstHalf << currentValue << secondHalf << std::endl;
				}
				needRerun = true;
			}
		} else {
			outStream << currentLine << std::endl;
		}
		std::getline(inStream, currentLine);
		if(currentLine.empty()){
			if(needRerun){
				inStream.close();
				outStream.close();
				inStream.open("__expfile.txt", std::ifstream::in);
				outStream.open("__tempfile.txt", std::ofstream::out);
				outStream << inStream.rdbuf();
				inStream.close();
				outStream.close();
				inStream.open("__tempfile.txt", std::ifstream::in);
				outStream.open("__expfile.txt", std::ofstream::out);
				needRerun = false;
				std::getline(inStream, currentLine);
			} else {
				inStream.close();
				outStream.close();
				break;
			}
		}
	}
	FILE* runfile = fopen("__expfile.txt", "r");
	c = fgetc(runfile);
	while(c != EOF){
		if(c == '\n') ++numberOfLines;
		c = fgetc(runfile);
	}
	rewind(runfile);
	runs = new mpf_class*[numberOfLines];
	maxOrders = new int[numberOfLines];
	for(int currentLine = 1; currentLine <= numberOfLines; ++currentLine){
		runs[currentLine-1] = new mpf_class[4];
		for(int i = 1; i <= 4; ++i){
			readStatus = ReadMPF(runs[currentLine-1][i-1], runfile);
			if(readStatus != 0) return -3;
		}
		maxOrders[currentLine-1] = ReadMaxOrder(runfile);
		if(maxOrders[currentLine-1] <= 0) return -3;
		fgetc(runfile);
	}
	std::remove("__tempfile.txt");
	std::remove("__expfile.txt");
	return numberOfLines;
}

int ReadMPF(mpf_class& output, FILE* runfile){
	ClearStructureChars(runfile);
	mpf_t temp;
	char inputBuffer[256]{0};
	int c = fgetc(runfile);
	int place = 0;
	while((c >= '0' && c <= '9') || c == '.'){
		inputBuffer[place] = (char)c;
		c = fgetc(runfile);
		++place;
	}
	if(place == 0) return -1;
	if(mpf_init_set_str(temp, inputBuffer, 10) == -1) return -1;
	mpf_set(output.get_mpf_t(), temp);
	mpf_clear(temp);
	ungetc(c, runfile);
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
	ClearStructureChars(runfile);
	int maxOrder = 0;
	int c = fgetc(runfile);
	while(c >= '0' && c <= '9'){
		maxOrder *= 10;
		maxOrder += c - '0';
		c = fgetc(runfile);
	}
	maxOrder -= maxOrder%2;
	ungetc(c, runfile);
	return maxOrder;
}

int RunCompare(mpf_class* run1, mpf_class* run2){
	if(run1[0] != run2[0]) return 0;	// runs are different
	if(run1[1] != run2[1]) return 0;
	if(run1[2] != run2[2]) return 0;
	if(run1[3] != run2[3]) return -2;	// runs differ only by hp
	return -1;							// runs are identical
}

void SetPowOverflow(unsigned short int maxOrder){
	powOverflow = new mpf_class[maxOrder/256+1];
	powOverflow[0] = 1;
	for(int i = 1; i <= maxOrder/256; ++i){
		powOverflow[i] = 16;
		mpf_pow_ui(powOverflow[i].get_mpf_t(), powOverflow[i].get_mpf_t(), 256*i);
	}
}

void DebugPrintRunVector(const mpf_class* runVector, const std::vector<mpf_class> hp, const unsigned short int maxOrder){
	for(unsigned int i = 1; i <= hp.size(); ++i){
		std::cout << "Pretend this is a run with c = " << to_string(runVector[0], 0) << ", hl = " << to_string(runVector[1], 0) << ", hh = " << to_string(runVector[2], 0) << ", hp = " << to_string(hp[i-1], 0) << ", maxOrder = " << maxOrder << "." << std::endl;
	}
}

void FindCoefficients(const mpf_class* runVector, const std::vector<mpf_class> hp, const unsigned short int maxOrder, const std::string outputName){
	int mnLocation[maxOrder]; 	/* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	int mnMultiplicity[maxOrder] = {0};	/* number of mn combinations giving i+1 = m*n */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int mTable[numberOfMN] = {0};
	unsigned short int nTable[numberOfMN] = {0};
	int mnLookup[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, numberOfMN, mnLocation, mnMultiplicity, maxOrder);
	
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	mpf_class bsq, invBsq, llsq, lhsq, temp1, temp2;
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);
	
	Cpqmn_t Cpqmn(&bsq, &invBsq, numberOfMN, maxOrder, mTable, nTable, mnLookup);
	Cpqmn.FillHpmn();
	
	auto time1 = Clock::now();
	Cpqmn.FillAmn();
	auto time2 = Clock::now();
	
	// combine Amn into Rmn
	time1 = Clock::now();
	Cpqmn.FillRmn(&llsq, &lhsq);
	time2 = Clock::now();
	
	// combine Rmn and hpmn into computation of H
	Hmn_t Hmn(&Cpqmn, numberOfMN, maxOrder, mnLocation, mnMultiplicity, mnLookup);
	time1 = Clock::now();
	Hmn.FillHmn();
	time2 = Clock::now();
	
	// corral H by q order and display coefficients
	mpf_class H[maxOrder/2+1];
	for(unsigned int i = 1; i <= hp.size(); ++i){
		H[0] = 1;
		FillH(H, &Hmn, &Cpqmn, hp[i-1], mnLocation, mnMultiplicity, maxOrder);
		if(outputName.empty()) DisplayH(H, runVector[0], runVector[1], runVector[2], hp[i-1], maxOrder);
		WriteH(H, runVector[0], runVector[1], runVector[2], hp[i-1], maxOrder, outputName);
	}
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

void ConvertInputs(mpf_class& bsq, mpf_class& invBsq, mpf_class& llsq, mpf_class& lhsq, const mpf_class& c, const mpf_class& hl, const mpf_class& hh, mpf_class& temp1, mpf_class& temp2){
	temp1 = c*c;
	temp2 = c*26;
	temp1 -= temp2;
	temp1 += 25;
	temp1 = sqrt(temp1);
	temp1 = c - temp1;
	temp1 -= 13;
	bsq = temp1/12;
	invBsq = 1/bsq;

	temp1 = c - 1;
	temp1 /= 24;
	llsq = hl - temp1;
	temp1 = c - 1;
	temp1 /= 24;
	lhsq = hh - temp1;
}

void FillH(mpf_class* H, const Hmn_t* Hmn, const Cpqmn_t* Cpqmn, const mpf_class hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	mpf_class temp1, temp2;
	for(int order = 2; order <= maxOrder; order+=2){
		for(int power = 2; power <= order; power+=2){
			for(int scanPos = mnLocation[power-1]; scanPos <= mnLocation[power-1] + mnMultiplicity[power-1] - 1; ++scanPos){
				temp1 = hp - Cpqmn->hpmn[scanPos-1];
				temp1 = Cpqmn->Rmn[scanPos-1]/temp1;
				temp1 *= Hmn->Hmn[(order-power)/2][scanPos-1];
				temp1 *= powOverflow[power/256];
				temp2 = 16;
				mpf_pow_ui(temp2.get_mpf_t(), temp2.get_mpf_t(), power%256);
				temp1 *= temp2;
				H[order/2] += temp1;
			}
		}
	}
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

std::string to_string(const mpf_class N, int digits){
	if(digits < 0) digits = -digits;	
	mp_exp_t dotPos;
	std::string output = N.get_str(dotPos, 10, digits);
	if(output.empty()) output.append("0");
	while(dotPos > (int)output.size()) output.append("0");
	if(digits == 0 && dotPos > 0 && dotPos < (int)output.size()){
		if(sgn(N) == -1) dotPos += 1;
		output.insert(dotPos, ".");
	}
	if(digits > 0 && dotPos > digits){
		while((int)output.size() < digits) output.append("0");
		if(sgn(N) == 1) output.insert(1, ".");		
		if(sgn(N) == -1) output.insert(2, ".");
		output.append("*10^");
		output.append(std::to_string(dotPos));
	}
	if(N > 0 && N < 1) output.insert(0, "0.");
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

void DisplayH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder){
	std::cout << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
	}
}

// An empty outputName means a single run; a filled one is a multirun, which prints fewer words.
void WriteH(const mpf_class* H, const mpf_class c, const mpf_class hl, const mpf_class hh, const mpf_class hp, const unsigned short int maxOrder, const std::string outputName){
	std::ofstream outputFile;
	std::string filename = outputName;
	if(outputName.empty()){
		filename = "virasoro_" + to_string(c, 3) + "_" + to_string(hl, 3) + "_" + to_string(hh, 3) + "_" + to_string(hp, 1) + "_" + std::to_string(maxOrder) + ".txt";
	}
	outputFile.open (filename, std::ios_base::app | std::ios_base::out);
	if(outputName.empty()){
		outputFile << "Given the parameters" << std::endl;
		outputFile << "c = " << to_string(c, 0) << ", h_L = " << to_string(hl, 0) << ", h_H = " << to_string(hh, 0) << ", h_p = " << to_string(hp, 0) << std::endl;
		outputFile << "the Virasoro block coefficients are as follows:" << std::endl;
		for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
			outputFile << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
		}
	} else {
		outputFile << "{" << to_string(c, 0) << "," << to_string(hl, 0) << "," << to_string(hh, 0) << "," << to_string(hp, 0) << "," << maxOrder << "}" << std::endl;
	}
	outputFile << "{1";
	for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
		outputFile << "," << to_string(H[orderBy2], 0);
	}
	outputFile << "}" << std::endl;
	outputFile.close();
}
