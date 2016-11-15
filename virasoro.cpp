#include "virasoro.h"

typedef std::chrono::high_resolution_clock Clock;

const int maxThreads = 8;			// Maximum number of simultaneous threads
const int precision = 512;			// Precision of mpf_class in bits
const mpf_class tolerance(1e-20);	// Smaller than this is taken to be 0 for comparisons
mpf_class* powOverflow;

int main(int argc, char** argv){
	mpf_set_default_prec(precision);
	auto programStart = Clock::now();
	std::string options = ParseOptions(argc, argv);
	int exitCode;
	switch(argc){
		default:	printf("Error: input either the name of a runfile or the five parameters c, hl, hh, hp, maxOrder\n");
					return EXIT_FAILURE;
		case 2:		if(options.find("m", 0) != std::string::npos){
						exitCode = RunFromFile(argv[1], options);
					} else {
						exitCode = RunFromFile(argv[1], options);
						ShowTime("Entire computation", programStart);
					}
					break;
		case 6:		if(options.find("m", 0) != std::string::npos){
						exitCode = RunFromTerminal(argv, options);
					} else {
						exitCode = RunFromTerminal(argv, options);
						ShowTime("Entire computation", programStart);
					}
					break;
	}

	return exitCode;
}

std::string ParseOptions(int &argc, char** &argv){
	std::string options = "";
	char** newArgv;
	int newArgc = argc;
	bool* realArg = new bool[argc]();
	for(int i = 1; i < argc; ++i){
		if(argv[i][0] != '-'){
			realArg[i] = true;
			continue;
		}
		if(strcmp(argv[i], "-m") == 0){
			options.append("m");
			--newArgc;
			realArg[i] = false;
		}
		if(strcmp(argv[i], "-c") == 0){
			options.append("c");
			--newArgc;
			realArg[i] = false;
		}
/*		if(strncmp(argv[i], "-p", 2) == 0){
			std::string precString = argv[i];
			precision = stoi(precString.substr(2), nullptr, 10);
			--newArgc;
			realArg[i] = false;
		}*/
	}
	newArgv = new char*[newArgc];
	newArgv[0] = argv[0];
	int pos = 1;
	for(int i = 1; i < argc; ++i){
		if(realArg[i]){
			newArgv[pos] = argv[i];
			++pos;
		}
	}
	argc = newArgc;
	argv = newArgv;
	delete[] newArgv;
	delete[] realArg;
	return options;
}

int RunFromFile(char* filename, const std::string options){
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;	
	mpf_class** rawRuns = (mpf_class**)malloc(sizeof(*rawRuns));
	int* rawMaxOrders = (int*)malloc(sizeof(*rawMaxOrders));
	int rawNumberOfRuns = ReadRunfile(filename, rawRuns, rawMaxOrders);
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
	bool* crunched = new bool[rawNumberOfRuns]();
	std::vector<mpf_class>* rawHp = new std::vector<mpf_class>[rawNumberOfRuns]();
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
	mpf_class** runs = new mpf_class*[numberOfRuns];
	std::vector<mpf_class>* hp = new std::vector<mpf_class>[numberOfRuns];
	int* maxOrders = new int[numberOfRuns];
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
	delete[] crunched;
	delete[] rawHp;
	std::free(rawRuns);
	std::free(rawMaxOrders);
	if(!wolframOutput){
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
	}
	int highestMax = 0;
	for(int i = 1; i <= numberOfRuns; ++i){
		if(maxOrders[i-1] > highestMax) highestMax = maxOrders[i-1];
	}
	SetPowOverflow(highestMax);
	auto runStart = Clock::now();
	std::string outputName;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
		for(int run = 1; run < numberOfRuns; ++run){
			FindCoefficients(runs[run-1], hp[run-1], maxOrders[run-1], outputName);
			std::cout << ",";
		}
		FindCoefficients(runs[numberOfRuns-1], hp[numberOfRuns-1], maxOrders[numberOfRuns-1], outputName);
		std::cout << "}";
	} else {
		if(consoleOutput){
			outputName = "__CONSOLE";
		} else {
			outputName = NameOutputFile(filename);
			std::remove(outputName.c_str());
		}
		for(int run = 1; run <= numberOfRuns; ++run){
			runStart = Clock::now();
//			DebugPrintRunVector(runs[run-1], hp[run-1], maxOrders[run-1]);
			printf("Beginning run %i of %i.\n", run, numberOfRuns);
			FindCoefficients(runs[run-1], hp[run-1], maxOrders[run-1], outputName);
			ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);
		}
	}
	delete[] runs;
	delete[] hp;
	delete[] maxOrders;
	return 0;
}

int RunFromTerminal(char** argv, const std::string options){
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;
	mpf_class runVector[4];
	for(int i = 1; i <= 4; ++i){
		runVector[i-1] = argv[i];
	}
	std::vector<mpf_class> hp;
	hp.push_back(runVector[3]);
	unsigned short int maxOrder = std::atoi(argv[5]);
	maxOrder -= (maxOrder % 2);
	SetPowOverflow(maxOrder);
//	DebugPrintRunVector(runVector, hp, maxOrder);
	std::string outputName;
	if(wolframOutput){
		outputName = "__MATHEMATICA";
		std::cout << "{";
	} else if(consoleOutput) {
		outputName = "__CONSOLE";
	} else {
		outputName = NameOutputFile(nullptr);
		std::remove(outputName.c_str());
	}
	FindCoefficients(runVector, hp, maxOrder, outputName);
	if(wolframOutput) std::cout << "}";
	return 0;
}

int ReadRunfile(char* filename, mpf_class** &runs, int* &maxOrders){
	int numberOfLines = 0;
	int readStatus = 0;
	int c;
	std::string workingString = ExpandRunFile(filename);
	const char* workingFilename = workingString.c_str();
	std::cout << "Expansion successful. Now working with filename " << workingFilename << "." << std::endl;
	if(strcmp(workingFilename, "__READFAIL") == 0) return -1;
	if(strcmp(workingFilename, "__NOBRACE") == 0) return -2;
	FILE* runfile = fopen(workingFilename,"r");
	c = fgetc(runfile);
	while(c != EOF){
		if(c == '\n') ++numberOfLines;
		c = fgetc(runfile);
	}
	rewind(runfile);
	std::cout << "There are " << numberOfLines << " lines." << std::endl;
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
	std::cout << "Lines processed successfully." << std::endl;
	std::remove(workingFilename);
	return numberOfLines;
}

std::string ExpandRunFile(char* filename){
	std::string expandedName, tempName;
	expandedName = filename;
	expandedName.append("__EXPANDED");
	ExpandBraces(filename);
	ExpandRelativeEqns(expandedName);
	ExpandBraces(expandedName);
	tempName = expandedName + "__EXPANDED";
	std::rename(tempName.c_str(), expandedName.c_str());
	return expandedName;
}

int ExpandBraces(std::string filename){
	std::string expandedName = filename+"__EXPANDED";
	std::string tempName = filename+"__TEMP";
	std::ifstream inStream;
	std::ofstream tempStream, outStream;
	inStream.open(filename, std::ifstream::in);
	outStream.open(expandedName, std::ofstream::out);
	if((inStream.rdstate() & std::ifstream::failbit) != 0 || (outStream.rdstate() & std::ofstream::failbit) != 0){
		return -1;
	}
	std::string currentLine, firstHalf, secondHalf, insideBraces;
	std::size_t leftPos, rightPos;
	mpf_class lowerBound, upperBound, increment;
	std::tuple<mpf_class, mpf_class, mpf_class> parsedBraces;
	std::getline(inStream, currentLine);
	bool needRerun = false;
	while(true){
		if((leftPos = currentLine.find("{")) != std::string::npos){
			firstHalf = currentLine.substr(0, leftPos);
			if((rightPos = currentLine.find("}", leftPos+1)) == std::string::npos){
				return -2;
			} else if(currentLine.find_first_not_of("0123456789-. ,;", leftPos+1) == rightPos) {
				secondHalf = currentLine.substr(rightPos+1, currentLine.length()-rightPos-1);
				insideBraces = currentLine.substr(leftPos+1, rightPos-leftPos-1);
				parsedBraces = ParseBraces(firstHalf, insideBraces);
				if(std::get<0>(parsedBraces) > std::get<1>(parsedBraces) || std::get<2>(parsedBraces) <= 0) return -2;
				for(mpf_class currentValue = std::get<0>(parsedBraces); currentValue <= std::get<1>(parsedBraces); currentValue += std::get<2>(parsedBraces)){
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
				inStream.open(expandedName, std::ifstream::in);
				outStream.open(tempName, std::ofstream::out);
				outStream << inStream.rdbuf();
				inStream.close();
				outStream.close();
				inStream.open(tempName, std::ifstream::in);
				outStream.open(expandedName, std::ofstream::out);
				needRerun = false;
				std::getline(inStream, currentLine);
			} else {
				inStream.close();
				outStream.close();
				break;
			}
		}
	}
	std::remove(tempName.c_str());
	return 0;
}

std::tuple<mpf_class, mpf_class, mpf_class> ParseBraces(std::string firstHalf, std::string insideBraces){
	mpf_class lowerBound, upperBound, increment;
	std::size_t numStart, numEnd;
	numStart = insideBraces.find_first_of("0123456789-.ch");
	numEnd = insideBraces.find_first_of(" ,;", numStart+1);
	if(insideBraces.find_first_of("ch") >= numEnd){
		lowerBound = insideBraces.substr(numStart, numEnd-numStart);
	} else {
		lowerBound = RelativeMPF(firstHalf, insideBraces.substr(numStart, numEnd-numStart));
	}
	numStart = insideBraces.find_first_of("0123456789-.ch", numEnd+1);
	numEnd = insideBraces.find_first_of(" ,;", numStart+1);
	if(insideBraces.find_first_of("ch") >= numEnd){
		upperBound = insideBraces.substr(numStart, numEnd-numStart);
	} else {
		upperBound = RelativeMPF(firstHalf, insideBraces.substr(numStart, numEnd-numStart));
	}
	numStart = insideBraces.find_first_of("0123456789-.ch", numEnd+1);
	numEnd = insideBraces.find_first_of(" ,;", numStart+1);
	if(insideBraces.find_first_of("ch") >= numEnd){
		increment = insideBraces.substr(numStart, numEnd-numStart);
	} else {
		increment = RelativeMPF(firstHalf, insideBraces.substr(numStart, numEnd-numStart));
	}
	return std::make_tuple(lowerBound, upperBound, increment);
}

void ExpandRelativeEqns(std::string filename){
	std::ifstream inStream;
	std::ofstream outStream;
	std::string currentLine, firstHalf, secondHalf;
	size_t leftPos, rightPos, splitPos;
	mpf_class value;
	std::string tempName = filename+"__TEMP";
	bool madeChange = false;
	inStream.open(filename);
	outStream.open(tempName);
	std::getline(inStream, currentLine);
	while(true){
		while(true){
			if((splitPos = currentLine.find_first_of("ch")) != std::string::npos){
				leftPos = currentLine.find_last_of(" ,;", splitPos-1);
				rightPos = currentLine.find_first_of(" ,;", splitPos+1);
				firstHalf = currentLine.substr(0, leftPos+1);
				secondHalf = currentLine.substr(rightPos, currentLine.length() - rightPos);
				value = RelativeMPF(firstHalf, currentLine.substr(leftPos+1, rightPos-leftPos-1));
				madeChange = true;
				outStream << firstHalf << value << secondHalf << std::endl;
			} else {
				outStream << currentLine << std::endl;
			}
			std::getline(inStream, currentLine);
			if(currentLine.empty()){
				inStream.close();
				outStream.close();
				break;
			}
		}
		if(madeChange){
			std::rename(tempName.c_str(), filename.c_str());
			inStream.open(filename);
			outStream.open(tempName);
			madeChange = false;
			std::getline(inStream, currentLine);
		} else {
			break;
		}
	}
	std::rename(tempName.c_str(), filename.c_str());
	return;
}

std::tuple<mpf_class, int> ParseRelativeEqn(std::string equation, std::string relTo){
	mpf_class modifier = 0;
	int type = -100;
	std::size_t hit;
	if((hit = equation.find(relTo)) != std::string::npos){
		if(hit == 0){
			if(equation[relTo.length()] == '+'){
				type = 0;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '-'){
				type = 1;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '*'){
				type = 3;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '/'){
				type = 4;
				modifier = equation.substr(relTo.size()+1);
			} else {				// assume it's just equality with no arithmetic
				type = 0;
				modifier = 0;
			}
		} else if(hit >= 2){
			if(equation[hit-1] == '+'){
				type = 0;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '-'){
				type = 2;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '*'){
				type = 3;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '/'){
				type = 5;
				modifier = equation.substr(0,hit-1);
			} else {				// assume it's just equality with no arithmetic
				type = 0;
				modifier = 0;
			}
		}
	}
	if(relTo == "c") type += 10;
	if(relTo == "hl") type += 20;
	if(relTo == "hh") type += 30;
	return std::make_tuple(modifier, type);
}

mpf_class RelativeMPF(std::string firstHalf, std::string equation){
	mpf_class output = 0;
	std::size_t baseStart, baseEnd;
	mpf_class baseMPF;
	std::tuple<mpf_class, int> parsedEqn;
	if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "c")) < 0){
		if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hl")) < 0){
			if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hh")) < 0){
				return 0;
			}
		}
	}
	mpf_class modifier = std::get<0>(parsedEqn);
	int type = std::get<1>(parsedEqn);
	switch(type){
		case 10:	// c + n
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF + modifier, 0);
					break;
		case 11:	// c - n
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF - modifier, 0);
					break;
		case 12:	// n - c
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier - baseMPF, 0);
					break;
		case 13:	// n*c
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF*modifier, 0);
					break;
		case 14:	// c/n
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF/modifier, 0);
					break;
		case 15:	// n/c
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier/baseMPF, 0);
					break;
		case 20:	// hl + n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF + modifier, 0);
					break;
		case 21:	// hl - n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF - modifier, 0);
					break;
		case 22:	// n - hl
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier - baseMPF, 0);
					break;
		case 23:	// n*hl
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF * modifier, 0);
					break;
		case 24:	// hl/n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF / modifier, 0);
					break;
		case 25:	// n/hl
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier / baseMPF, 0);
					break;
		case 30:	// hh + n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier + baseMPF, 0);
					break;
		case 31:	// hh - n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF - modifier, 0);
					break;
		case 32:	// n - hh
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier - baseMPF, 0);
					break;
		case 33:	// n*hh
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier * baseMPF, 0);
					break;
		case 34:	// hh/n
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(baseMPF/modifier, 0);
					break;
		case 35:	// n/hh
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					firstHalf.erase(0, firstHalf.find_first_of(" ,;")+1);
					baseStart = 0;
					baseEnd = firstHalf.find_first_of(" ,;");
					baseMPF = firstHalf.substr(baseStart, baseEnd - baseStart);
					output = to_string(modifier/baseMPF, 0);
					break;
	}
	return output;
}

int RunCompare(mpf_class* run1, mpf_class* run2){
	if(run1[0] != run2[0]) return 0;
	if(run1[1] != run2[1]) return 0;	// runs are different
	if(run1[2] != run2[2]) return 0;
	if(run1[3] != run2[3]) return -2;	// runs differ only by hp
	return -1;							// runs are identical
}

int ReadMPF(mpf_class& output, FILE* runfile){
	ClearStructureChars(runfile);
	mpf_t temp;
	char inputBuffer[256]{0};
	int c = fgetc(runfile);
	int place = 0;
	while((c >= '0' && c <= '9') || c == '.' || c == '-'){
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

void FindCoefficients(const mpf_class* runVector, const std::vector<mpf_class> hp, unsigned short int maxOrder, const std::string outputName){
	// construct b^2 and 1/b^2 from c and lambda_l and lambda_h from h_l and h_h
	mpf_class bsq, invBsq, llsq, lhsq, temp1, temp2;
	ConvertInputs(bsq, invBsq, llsq, lhsq, runVector[0], runVector[1], runVector[2], temp1, temp2);
	CheckForDivergences(&invBsq, maxOrder);
	if(maxOrder <= 2) return;

	int* mnLocation = new int[maxOrder]; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	int* mnMultiplicity = new int[maxOrder]();	/* number of mn combinations giving i+1 = m*n */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int* mTable = new unsigned short int[numberOfMN]();
	unsigned short int* nTable = new unsigned short int[numberOfMN]();
	int* mnLookup = new int[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, mnLocation, mnMultiplicity, maxOrder);
	
	Cpqmn_t Cpqmn(&bsq, &invBsq, numberOfMN, maxOrder, mTable, nTable, mnLookup);
	Cpqmn.FillHpmn();
	
	auto time1 = Clock::now();
	Cpqmn.FillRmn(&llsq, &lhsq);
	auto time2 = Clock::now();
	
	// combine Rmn and hpmn into computation of H
	Hmn_t Hmn(&Cpqmn, numberOfMN, maxOrder, mnLocation, mnMultiplicity, mnLookup);
	time1 = Clock::now();
	Hmn.FillHmn();
	time2 = Clock::now();
	
	// corral H by q order and display coefficients
	mpf_class* H = new mpf_class[maxOrder/2+1];
	for(unsigned int i = 1; i <= hp.size(); ++i){
		H[0] = 1;
		FillH(H, &Hmn, &Cpqmn, hp[i-1], mnLocation, mnMultiplicity, maxOrder);
		if(outputName.empty() || outputName == "__CONSOLE") DisplayH(H, runVector[0], runVector[1], runVector[2], hp[i-1], maxOrder);
		if(outputName != "__CONSOLE") WriteH(H, runVector[0], runVector[1], runVector[2], hp[i-1], maxOrder, outputName);
	}
	delete[] mnLocation;
	delete[] mnMultiplicity;
	delete[] mTable;
	delete[] nTable;
	delete[] mnLookup;
	delete[] H;
	return;
}

void CheckForDivergences(const mpf_class* bsq, unsigned short int &maxOrder){
	mpf_class temp;
	int lowest = maxOrder;
	for(int m = 1; m <= maxOrder-2; ++m){
		for(int n = 1+(m%2); m*n <= maxOrder-2; n+=1+(m%2)){
			for(int p = 1; p <= maxOrder-m*n; ++p){
				for(int q = 1+(p%2); p*q <= maxOrder-m*n; q+=1+(p%2)){
					temp = (*bsq)*(n-q) - (m+p);
					if(abs(temp) < tolerance){
						if(m*n < lowest) lowest = m*n;
//						if(p*q < lowest) lowest = p*q;
						continue;
					}
					temp = (*bsq)*(n+q) - (m-p);
					if(abs(temp) < tolerance){
						if(m*n < lowest) lowest = m*n;
//						if(p*q < lowest) lowest = p*q;
						continue;
					}
				}
			}
		}
	}
	if(*bsq == 2) lowest = 2; // this is a shitty workaround until I can fix the 1/0 in Amn
	lowest = lowest - (lowest%2);
	if(lowest < maxOrder){
		maxOrder = lowest;
		if(lowest > 2) printf("Stopping this run at order %i because the coefficients diverge above this.\n",lowest);
		if(lowest <= 2) printf("Skipping this run because the coefficients diverge immediately.\n");
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
	if(digits > 0 && dotPos < digits){
		while(dotPos > (int)output.size()) output.append("0");
		if((N > 1 || N < 0) && dotPos < (int)output.size()) output.insert(abs(dotPos), ".");
	}
	if(digits == 0 && dotPos > 0){
		while(dotPos > (int)output.size()) output.append("0");
		if(dotPos < (int)output.size()){
			if(sgn(N) == -1) dotPos += 1;
			output.insert(dotPos, ".");
		}
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
	if(outputName == "__MATHEMATICA"){
		std::cout << "{" << to_string(c, 0) << "," << to_string(hl, 0) << "," << to_string(hh, 0) << "," << to_string(hp, 0) << "," << maxOrder << "},{1";
		for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
			std::cout << "," << to_string(H[orderBy2], 0);
		}
		std::cout << "}";
		return;
	}
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
