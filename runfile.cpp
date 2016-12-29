#include "runfile.h"

Runfile_c::Runfile_c(){}
Runfile_c::Runfile_c(const Runfile_c& other): filename(other.filename), lines(other.lines), maxThreads(other.maxThreads), precision(other.precision), tolerance(other.tolerance), showProgressBar(other.showProgressBar){}
Runfile_c::Runfile_c(Runfile_c&& other): Runfile_c(){	swap(*this, other);	}
Runfile_c::Runfile_c(const char* filename): filename(filename)
{
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	std::string currentLine;
	while(true){
		std::getline(inStream, currentLine);
		if(currentLine.empty()) break;
		lines.push_back(currentLine);
	}
}
Runfile_c::Runfile_c(const std::string filename): filename(filename)
{
	std::ifstream inStream;
	inStream.open(filename, std::ifstream::in);
	std::string currentLine;
	while(true){
		std::getline(inStream, currentLine);
		if(currentLine.empty()) break;
		lines.push_back(currentLine);
	}
}
Runfile_c::Runfile_c(const std::vector<std::string> line): filename("command_line"){
	std::string combinedLine = "";
	for(unsigned int i = 1; i <= line.size(); ++i){
		combinedLine += line[i-1];
		combinedLine += " ";
	}
	combinedLine.erase(combinedLine.length()-1);
	lines.push_back(combinedLine);
}
Runfile_c::~Runfile_c(){
}
void Runfile_c::swap(Runfile_c& first, Runfile_c& second){
	std::swap(first.filename, second.filename);
	std::swap(first.lines, second.lines);
	std::swap(first.maxThreads, second.maxThreads);
	std::swap(first.precision, second.precision);
	std::swap(first.tolerance, second.tolerance);
	std::swap(first.showProgressBar, second.showProgressBar);
	return;
}

Runfile_c& Runfile_c::operator=(Runfile_c&& v){
	swap(*this, v);
	return *this;
}
Runfile_c& Runfile_c::operator=(Runfile_c v){
	swap(*this, v);
	return *this;
}
Runfile_c& Runfile_c::operator=(const char* &filename){
	Runfile_c newFile(filename);
	swap(*this, newFile);
	return *this;
}
Runfile_c& Runfile_c::operator=(const std::string& filename){
	Runfile_c newFile(filename);
	swap(*this, newFile);
	return *this;
}
Runfile_c& Runfile_c::operator=(const std::vector<std::string> line){
	Runfile_c newFile(line);
	swap(*this, newFile);
	return *this;
}

void Runfile_c::SetMaxThreads(int newMax){
	maxThreads = newMax;
	return;
}
void Runfile_c::SetPrecision(int newPrec){
	precision = newPrec;
	return;
}
void Runfile_c::SetTolerance(mpf_class newTolerance){
	tolerance = newTolerance;
	return;
}
void Runfile_c::SetProgressBar(bool newProgressBar){
	showProgressBar = newProgressBar;
	return;
}
int Runfile_c::NumberOfRuns(){
	return (int)lines.size();	// this should be runs.size() after expanding but whatever
}

int Runfile_c::ReadRunfile(){
	if(lines.empty()){
		perror("Error: if you give one argument it must be the name of a runfile.\n");
		return -1;
	}
	int errorCode = Expand();
	if(errorCode == -2){
		perror("Error: a curly brace '{' was found indicating a batch run, but it was not followed by a valid macro.\n");
		return -2;
	}
	size_t leftPos, rightPos;
	std::vector<mpfc_class> currentRun;
	int currentMO;
	for(unsigned int i = 1; i <= lines.size(); ++i){
//		std::cout << "About to parse the following line:" << std::endl;
//		std::cout << lines[i-1] << std::endl;
		leftPos = 0;
		rightPos = 0;
		for(int j = 1; j <= 4; ++j){
			leftPos = lines[i-1].find_first_of("0123456789.-(", rightPos);
			if(lines[i-1][leftPos] == '('){
				rightPos = lines[i-1].find(")", leftPos) + 1;
			} else {
				rightPos = lines[i-1].find_first_not_of("0123456789.-", leftPos);
			}
//			std::cout << "Going to emplace with this:" << lines[i-1].substr(leftPos, rightPos-leftPos) << std::endl;
			currentRun.emplace_back(lines[i-1].substr(leftPos, rightPos-leftPos));
/*			if(emplacing doesn't work){
				perror("Error: expected a number in the runfile but failed to read one.\n");
				return -3;
			}*/
		}
		runs.push_back(currentRun);
		leftPos = lines[i-1].find_first_of("0123456789.-", rightPos);
		rightPos = lines[i-1].find_first_not_of("0123456789.-", leftPos);
		currentMO = std::stoi(lines[i-1].substr(leftPos, rightPos-leftPos));
		maxOrders.push_back(currentMO);
		currentRun.clear();
	}
	if(runs.empty()){
		perror("Error: zero valid runs detected in given runfile.\n");
		return 0;
	}
	// Check runs for duplicates, compress runs differing only by hp
	bool* crunched = new bool[runs.size()]();
	for(unsigned int i = 1; i <= runs.size(); ++i){
		for(unsigned int j = i+1; j <= runs.size(); ++j){
			if(crunched[j-1]) continue;
			switch(RunCompare(runs[i-1], runs[j-1])){
				case 0:		crunched[j-1] = false;	// runs different, do nothing
							break;				
				case -1:	if(maxOrders[j-1] > maxOrders[i-1]) maxOrders[i-1] = maxOrders[j-1];
							crunched[j-1] = true;	// runs identical, crunch them together
							break;
				case -2:	for(unsigned int k = 4; k <= runs[j-1].size(); ++k) runs[i-1].push_back(runs[j-1][k-1]);
							if(maxOrders[j-1] > maxOrders[i-1]) maxOrders[i-1] = maxOrders[j-1];
							crunched[j-1] = true;	// runs differ by hp, make them fast
							break;
			}
		}
	}
	for(unsigned int i = runs.size(); i >= 1; --i){
		if(crunched[i-1]){
			runs.erase(runs.begin()+i-1);
			maxOrders.erase(maxOrders.begin()+i-1);
		}
	}
	delete[] crunched;
	return lines.size();
}

int Runfile_c::Expand(){
	for(int param = 1; param <= 4; ++param){
		ExpandRelativeEqns(param-1);
		ExpandBraces(param-1);
	}
	return 0; 
}

// param=0 is c/b/b^2; param=1 is hl; param=2 is hh; param=3 is hp
int Runfile_c::ExpandBraces(const int param){
	std::vector<std::string> newLines;
	std::tuple<mpfc_class, mpfc_class, mpfc_class> parsedBraces;
	mpfc_class currentValue;
	std::string currentParam, insideBraces, firstHalf, secondHalf;
	size_t paramLeftPos, paramRightPos, leftPos, rightPos;
	std::tuple<size_t, size_t> paramLocation;
	mpfc_class value;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		paramLocation = FindNthParameter(lines[i-1], param);
		paramLeftPos = std::get<0>(paramLocation);
		paramRightPos = std::get<1>(paramLocation);
		currentParam = lines[i-1].substr(paramLeftPos, paramRightPos-paramLeftPos+1);
//		std::cout << "Checking " << currentParam << " for braces." << std::endl;
		if((leftPos = currentParam.find("{")) != std::string::npos){
			firstHalf = lines[i-1].substr(0, paramLeftPos+leftPos);
			if((rightPos = currentParam.find("}", leftPos+1)) == std::string::npos){
//				std::cout << "Found unpaired braces in " << currentParam << "." << std::endl;
				return -2;
			} else if(currentParam.find_first_not_of("0123456789-+()e. ,;", leftPos+1) == rightPos) {
//				std::cout << "Found paired braces in " << currentParam << ", parsing." << std::endl;
				secondHalf = lines[i-1].substr(paramLeftPos+rightPos+1);
				insideBraces = currentParam.substr(leftPos+1, rightPos-leftPos-1);
				parsedBraces = ParseBraces(insideBraces);
//				if(std::get<0>(parsedBraces).realPart() => std::get<1>(parsedBraces).realPart() || std::get<2>(parsedBraces).realPart() <= 0) return -2;
//				std::cout << "Parsed " << currentParam << " into the following lines:" << std::endl;
				for(currentValue = std::get<0>(parsedBraces); currentValue.realPart() <= std::get<1>(parsedBraces).realPart(); currentValue += std::get<2>(parsedBraces)){
					newLines.push_back(firstHalf + to_string(currentValue, -1) + secondHalf);
//					std::cout << firstHalf + to_string(currentValue, -1) + secondHalf << std::endl;
				}
			} else {
//				std::cout << "Got confused by " << currentParam << ", skipping." << std::endl;
				newLines.push_back(lines[i-1]);
			}
		} else {
			newLines.push_back(lines[i-1]);
		}
	}
	lines = newLines;
	newLines.clear();
	return 0;
}

std::tuple<mpfc_class, mpfc_class, mpfc_class> Runfile_c::ParseBraces(std::string insideBraces){
	std::size_t numStart, numEnd;

	numStart = insideBraces.find_first_of("(0123456789-+.ch");
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	mpfc_class lowerBound(insideBraces.substr(numStart, numEnd-numStart+1));
//	std::cout << "Lower bound is " << to_string(lowerBound, 4) << " parsed from " << insideBraces.substr(numStart, numEnd-numStart+1) << std::endl;

	numStart = insideBraces.find_first_of("(0123456789-+.ch", numEnd+2);
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	mpfc_class upperBound(insideBraces.substr(numStart, numEnd-numStart+1));
//	std::cout << "Upper bound is " << to_string(upperBound, 4) << " parsed from " << insideBraces.substr(numStart, numEnd-numStart+1) << std::endl;

	numStart = insideBraces.find_first_of("(0123456789-+.ch", numEnd+2);
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	mpfc_class increment(insideBraces.substr(numStart, numEnd-numStart+1));
//	std::cout << "Increment is " << to_string(increment, 4) << " parsed from " << insideBraces.substr(numStart, numEnd-numStart+1).c_str() << std::endl;

	return std::make_tuple(lowerBound, upperBound, increment);
}

std::tuple<size_t, size_t> Runfile_c::FindNthParameter(const std::string line, const int param){
//	std::cout << "Finding parameter #" << param << " from the following line: " << line << std::endl;
	size_t splitPos;
	size_t leftPos = 0;
	size_t rightPos = 0;
	int paramsFound = 0;
	std::vector<size_t> paramLocations;
	paramLocations.push_back(0);
	do{
		splitPos = line.find_first_of(" ", rightPos+2);
		leftPos = line.find_last_of("(", splitPos-1);
		rightPos = line.find_last_of(")", splitPos-1);
		if(leftPos == std::string::npos || (rightPos!=std::string::npos && leftPos < rightPos)){
			rightPos = splitPos-1;
			++paramsFound;
			paramLocations.push_back(rightPos+2);
		} else {
			rightPos = splitPos-1;
		}
	}while(paramsFound <= param);
//	std::cout << "It's " << line.substr(paramLocations[param], paramLocations[param+1]-2-paramLocations[param]+1) << " between positions " << paramLocations[param] << " and " << paramLocations[param+1] << std::endl;
	return std::make_tuple(paramLocations[param], paramLocations[param+1]-2);
}

// param=0 is c/b/b^2; param=1 is hl; param=2 is hh; param=3 is hp
int Runfile_c::ExpandRelativeEqns(const int param){
/*	std::cout << "Parsing parameter " << param << " from these lines:" << std::endl;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		std::cout << lines[i-1] << std::endl;
	}*/
	int changesMade = 0;
	std::vector<std::string> newLines;
	std::string currentParam, newParam, firstHalf, secondHalf;
	size_t paramLeftPos, paramRightPos, leftPos, rightPos, splitPos;
	std::tuple<size_t,size_t> paramLocation;
	mpfc_class value;
	bool madeChange = true;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		paramLocation = FindNthParameter(lines[i-1], param);
		paramLeftPos = std::get<0>(paramLocation);
		paramRightPos = std::get<1>(paramLocation);
		currentParam = lines[i-1].substr(paramLeftPos, paramRightPos-paramLeftPos+1);
//		std::cout << "Checking " << currentParam << " for relative MPFs." << std::endl;
		do{
			madeChange = false;
			if((splitPos = currentParam.find_first_of("ch")) != std::string::npos){
				leftPos = currentParam.find_last_of(" ,;{", splitPos-1) + 1;
				rightPos = currentParam.find_first_of(" ,;}", splitPos+1) - 1;
				if(rightPos >= currentParam.length()) rightPos = currentParam.length()-1;
				firstHalf = currentParam.substr(0, leftPos);
				secondHalf = currentParam.substr(rightPos+1);
				value = RelativeMPF(lines[i-1].substr(0, paramLeftPos), currentParam.substr(leftPos, rightPos-leftPos+1));
				madeChange = true;
				++changesMade;
				currentParam = firstHalf + to_string(value, -1) + secondHalf;
//				std::cout << currentParam << "." << std::endl;
			}
		}while(madeChange);
		newLines.push_back(lines[i-1].substr(0, paramLeftPos) + currentParam + lines[i-1].substr(paramRightPos+1));
	}
	lines = newLines;
	newLines.clear();
/*	std::cout << "Ended up with these lines:" << std::endl;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		std::cout << lines[i-1] << std::endl;
	}*/
	return changesMade;
}

std::tuple<mpfc_class, int> Runfile_c::ParseRelativeEqn(std::string equation, std::string relTo){
	mpfc_class modifier = 0;
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

mpfc_class Runfile_c::RelativeMPF(std::string firstHalf, std::string equation){
	mpfc_class output = 0;
	mpfc_class baseMPF;
	std::tuple<mpfc_class, int> parsedEqn;
	if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "c")) < 0){
		if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hl")) < 0){
			if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hh")) < 0){
				return 0;
			}
		}
	}
	mpfc_class modifier = std::get<0>(parsedEqn);
	int type = std::get<1>(parsedEqn);
	switch(type){
		case 10:	// c + n
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = baseMPF + modifier;
					break;
		case 11:	// c - n
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = baseMPF - modifier;
					break;
		case 12:	// n - c
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = modifier - baseMPF;
					break;
		case 13:	// n*c
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = baseMPF*modifier;
					break;
		case 14:	// c/n
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = baseMPF/modifier;
					break;
		case 15:	// n/c
					baseMPF = FindBaseNumber(firstHalf, 0);
					output = modifier/baseMPF;
					break;
		case 20:	// hl + n
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = baseMPF + modifier;
					break;
		case 21:	// hl - n
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = baseMPF - modifier;
					break;
		case 22:	// n - hl
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = modifier - baseMPF;
					break;
		case 23:	// n*hl
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = baseMPF * modifier;
					break;
		case 24:	// hl/n
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = baseMPF / modifier;
					break;
		case 25:	// n/hl
					baseMPF = FindBaseNumber(firstHalf, 1);
					output = modifier / baseMPF;
					break;
		case 30:	// hh + n
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = modifier + baseMPF;
					break;
		case 31:	// hh - n
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = baseMPF - modifier;
					break;
		case 32:	// n - hh
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = modifier - baseMPF;
					break;
		case 33:	// n*hh
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = modifier * baseMPF;
					break;
		case 34:	// hh/n
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = baseMPF/modifier;
					break;
		case 35:	// n/hh
					baseMPF = FindBaseNumber(firstHalf, 2);
					output = modifier/baseMPF;
					break;
	}
	return output;
}

std::string Runfile_c::FindBaseNumber(std::string sourceString, const int paramNumber){
	std::size_t baseStart, baseEnd;
	baseStart = 0;
	for(int i = 1; i <= paramNumber; ++i){
		if(sourceString[baseStart] == '('){
			sourceString.erase(0, sourceString.find(")")+2);
		} else {
			sourceString.erase(0, sourceString.find_first_of(" ,;")+1);
		}
	}
	if(sourceString[baseStart] == '('){
		baseEnd = sourceString.find(")")+1;
	} else {
		baseEnd = sourceString.find_first_of(" ,;");
	}
	return sourceString.substr(baseStart, baseEnd - baseStart);
}

int Runfile_c::RunCompare(std::vector<mpfc_class> run1, std::vector<mpfc_class> run2){
/*	std::cout << "Comparing these two runs:" << std::endl;
	for(unsigned int i = 1; i <= run1.size(); ++i) std::cout << run1[i-1] << " ";
	std::cout << std::endl;
	for(unsigned int i = 1; i <= run2.size(); ++i) std::cout << run2[i-1] << " ";
	std::cout << std::endl;*/
	if((run1[0] != run2[0]) || (run1[1] != run2[1]) || (run1[2] != run2[2])) return 0;	// runs are different
	for(unsigned int i = 4; i <= run1.size(); ++i){
		for(unsigned int j = 4; j <= run2.size(); ++j){
			if(run1[i-1] != run2[j-1]) return -2;	// runs differ only by hp
		}
	}
	return -1;							// runs are identical
}

std::string Runfile_c::NameOutputFile(){
	std::string outputname = filename;
	if(lines.size() == 1){
		outputname = "virasoro_" + to_string(runs[0][0], 3) + "_" + to_string(runs[0][1], 3) + "_" + to_string(runs[0][2], 3) + "_" + to_string(runs[0][3], 1) + "_" + std::to_string(maxOrders[0]) + ".txt";
	} else {
		std::size_t delPos = outputname.find(".txt");
		if(delPos != std::string::npos) outputname.erase(delPos, 4);
		outputname.append("_results.txt");
	}
	return outputname;
}

int Runfile_c::Execute(std::string options){
	auto programStart = Clock::now();
	const bool wolframOutput = options.find("m", 0) != std::string::npos;
	const bool consoleOutput = options.find("c", 0) != std::string::npos;	
	const bool wstp = options.find("w", 0) != std::string::npos;
	int bGiven = 0;
	if(options.find("b", 0) != std::string::npos) bGiven = 1;
	if(options.find("bb", 0) != std::string::npos) bGiven = 2;
	if(ReadRunfile() <= 0) return -1;
	if(!wstp && !wolframOutput){
		for(unsigned int i = 1; i <= runs.size(); ++i){
			std::cout << "Run " << i << ": ";
			for(int j = 1; j <= 3; ++j){
				std::cout << to_string(runs[i-1][j-1], 4) << " ";
			}
			if(runs[i-1].size() > 4) std::cout << "{";
			for(unsigned int j = 4; j <= runs[i-1].size(); ++j){
				std::cout << to_string(runs[i-1][j-1], 4) << ",";
			}
			if(runs[i-1].size() > 4){
				std::cout << "\b} ";
			} else {
				std::cout << "\b ";
			}
			std::cout << maxOrders[i-1];
			std::cout << std::endl;
		}
		if(!consoleOutput) std::cout << "Output will be saved to " << NameOutputFile() << ". If it exists, it will be overwritten." << std::endl;
	}
	int highestMax = 0;
	for(unsigned int i = 1; i <= runs.size(); ++i){
		if(maxOrders[i-1] > highestMax) highestMax = maxOrders[i-1];
	}
	auto runStart = Clock::now();
	std::string outputName;
	std::vector<mpf_class> realRunVector;
	bool allReal;
	if(wstp || wolframOutput){
		showProgressBar = false;
		if(wolframOutput){
			outputName = "__MATHEMATICA";
			std::cout << "{";
		} else {
			outputName = "__WSTP";
		}
		for(unsigned int run = 1; run <= runs.size(); ++run){
			allReal = true;
			for(unsigned int i = 1; i <= runs[run-1].size(); ++i){
				if(!runs[run-1][i-1].isReal()){
					allReal = false;
					break;
				}
			}
			if(bGiven == 0 && runs[run-1][0].realPart() < 25 && runs[run-1][0].realPart() > 1) allReal = false;
			if(allReal){
				for(unsigned int i = 1; i <= runs[run-1].size(); ++i) realRunVector.push_back(runs[run-1][i-1].realPart());
				FindCoefficients<mpf_class>(realRunVector, maxOrders[run-1], outputName, bGiven);
				realRunVector.clear();
			} else {
				FindCoefficients<mpfc_class>(runs[run-1], maxOrders[run-1], outputName, bGiven);
			}
			if(run < runs.size()) std::cout << ",";
		}
		if(wolframOutput) std::cout << "}";
	} else {
		if(consoleOutput){
			outputName = "__CONSOLE";
		} else {
			outputName = NameOutputFile();
			std::remove(outputName.c_str());
		}
		for(unsigned int run = 1; run <= runs.size(); ++run){
			runStart = Clock::now();
//			DebugPrintRunVector(runs[run-1], hp[run-1], maxOrders[run-1]);
			std::cout << "Beginning run " << run << " of " << runs.size() << "." << std::endl;
			allReal = true;
			for(unsigned int i = 1; i <= runs[run-1].size(); ++i){
				if(!runs[run-1][i-1].isReal()){
					allReal = false;
					break;
				}
			}
			if(bGiven == 0 && runs[run-1][0].realPart() < 25 && runs[run-1][0].realPart() > 1) allReal = false;
			if(allReal){
				for(unsigned int i = 1; i <= runs[run-1].size(); ++i) realRunVector.push_back(runs[run-1][i-1].realPart());
				FindCoefficients<mpf_class>(realRunVector, maxOrders[run-1], outputName, bGiven);
				realRunVector.clear();
			} else {
				FindCoefficients<mpfc_class>(runs[run-1], maxOrders[run-1], outputName, bGiven);
			}
			if(runs.size() > 1) ShowTime(std::string("Computing run ").append(std::to_string(run)), runStart);
		}
	}
	if(!wstp && !wolframOutput) ShowTime("Entire computation", programStart);
	return 0;
}

int Runfile_c::EnumerateMN (int* mnLocation, int* mnMultiplicity, const unsigned short int maxOrder){
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

void Runfile_c::FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
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

void Runfile_c::ShowTime(std::string computationName, std::chrono::time_point<std::chrono::high_resolution_clock> timeStart){
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
	mp_exp_t expo, dotPos;
	std::string output = N.get_str(expo, 10, digits);
	dotPos = expo;
	--expo;
	if(dotPos >= 0 && sgn(N) == -1) ++dotPos;
	if(dotPos < 0 && sgn(N) == -1) --dotPos;
	double Nd = N.get_d();
	if(output.empty()) return "0";
	if(digits > 0 && abs(expo) <= digits-1){			// number small enough, just write it
		while(dotPos > (int)output.size()) output.append("0");
		if(dotPos < (int)output.size()){
			if(Nd >= 1 || Nd <= -1){
				output.insert(dotPos, ".");
			} else if(Nd < 1 && Nd > 0) {
				output.insert(0, std::max(-(int)dotPos,0), '0');
				output.insert(0, "0.");
			} else if(Nd < 0 && Nd > -1){
				output.insert(1, std::max(-(int)dotPos,0), '0');
				output.insert(1, "0.");
			}
		}
	}
	if(digits > 0 && abs(expo) > digits-1){			// number too big, use a*10^b
		while((int)output.size() < digits) output.append("0");
		if(sgn(N) == 1) output.insert(1, ".");		
		if(sgn(N) == -1){
			output.append("0");
			output.insert(2, ".");
		}
		output.append("*10^");
		output.append(std::to_string(expo));
	}
	if(digits == 0){								// entire number has been requested
		if(Nd >= 1 || Nd <= -1){
			while(dotPos > (int)output.size()) output.append("0");
			if(dotPos < (int)output.size()) output.insert(dotPos, ".");
		} else if(Nd > 0) {
			output.insert(0, std::max(-(int)dotPos, 0), '0');
			output.insert(0, "0.");
		} else if(Nd > -1) {
			output.insert(1, std::max(-(int)dotPos, 0), '0');
			output.insert(1, "0.");
		}
	}
	return output;
}
