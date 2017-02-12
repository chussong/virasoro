#include "runfile.h"
#include "virasoro.h"

namespace virasoro {

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

int Runfile_c::NumberOfRuns(){
	int count = 0;
	for(unsigned int i = 1; i <= runs.size(); ++i){
		count += runs[i-1].size()-3;
	}
	return count;
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
	std::vector<std::complex<mpfr::mpreal>> currentRun;
	int currentMO;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		leftPos = 0;
		rightPos = 0;
		for(int j = 1; j <= 4; ++j){
			leftPos = lines[i-1].find_first_of("0123456789.-(", rightPos);
			if(lines[i-1][leftPos] == '('){
				rightPos = lines[i-1].find(")", leftPos) + 1;
			} else {
				rightPos = lines[i-1].find_first_not_of("0123456789.-", leftPos);
			}
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
	std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> parsedBraces;
	std::complex<mpfr::mpreal> currentValue;
	std::string currentParam, insideBraces, firstHalf, secondHalf;
	size_t paramLeftPos, paramRightPos, leftPos, rightPos;
	std::tuple<size_t, size_t> paramLocation;
	std::complex<mpfr::mpreal> value;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		paramLocation = FindNthParameter(lines[i-1], param);
		paramLeftPos = std::get<0>(paramLocation);
		paramRightPos = std::get<1>(paramLocation);
		currentParam = lines[i-1].substr(paramLeftPos, paramRightPos-paramLeftPos+1);
		if((leftPos = currentParam.find("{")) != std::string::npos){
			firstHalf = lines[i-1].substr(0, paramLeftPos+leftPos);
			if((rightPos = currentParam.find("}", leftPos+1)) == std::string::npos){
				return -2;
			} else if(currentParam.find_first_not_of("0123456789-+()e. ,;", leftPos+1) == rightPos) {
				secondHalf = lines[i-1].substr(paramLeftPos+rightPos+1);
				insideBraces = currentParam.substr(leftPos+1, rightPos-leftPos-1);
				parsedBraces = ParseBraces(insideBraces);
				for(currentValue = std::get<0>(parsedBraces); currentValue.real() <= std::get<1>(parsedBraces).real(); currentValue += std::get<2>(parsedBraces)){
					newLines.push_back(firstHalf + to_string(currentValue, -1) + secondHalf);
				}
			} else {
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

std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> Runfile_c::ParseBraces(std::string insideBraces){
	std::size_t numStart, numEnd;

	numStart = insideBraces.find_first_of("(0123456789-+.ch");
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	std::complex<mpfr::mpreal> lowerBound(insideBraces.substr(numStart, numEnd-numStart+1));

	numStart = insideBraces.find_first_of("(0123456789-+.ch", numEnd+2);
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	std::complex<mpfr::mpreal> upperBound(insideBraces.substr(numStart, numEnd-numStart+1));

	numStart = insideBraces.find_first_of("(0123456789-+.ch", numEnd+2);
	if(insideBraces[numStart] == '('){
		numEnd = insideBraces.find_first_of(")", numStart+1);
	} else {
		numEnd = insideBraces.find_first_of(" ,;", numStart+1)-1;
	}
	std::complex<mpfr::mpreal> increment(insideBraces.substr(numStart, numEnd-numStart+1));

	return std::make_tuple(lowerBound, upperBound, increment);
}

std::tuple<size_t, size_t> Runfile_c::FindNthParameter(const std::string line, const int param){
	size_t splitPos;
	size_t leftPos;
	size_t rightPos = -1;
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
	std::complex<mpfr::mpreal> value;
	bool madeChange = true;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		paramLocation = FindNthParameter(lines[i-1], param);
		paramLeftPos = std::get<0>(paramLocation);
		paramRightPos = std::get<1>(paramLocation);
		currentParam = lines[i-1].substr(paramLeftPos, paramRightPos-paramLeftPos+1);
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

std::tuple<std::complex<mpfr::mpreal>, int> Runfile_c::ParseRelativeEqn(std::string equation, std::string relTo){
	std::complex<mpfr::mpreal> modifier(0);
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

std::complex<mpfr::mpreal> Runfile_c::RelativeMPF(std::string firstHalf, std::string equation){
	std::complex<mpfr::mpreal> output(0);
	std::complex<mpfr::mpreal> baseMPF;
	std::tuple<std::complex<mpfr::mpreal>, int> parsedEqn;
	if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "c")) < 0){
		if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hl")) < 0){
			if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hh")) < 0){
				return std::complex<mpfr::mpreal>(0);
			}
		}
	}
	std::complex<mpfr::mpreal> modifier = std::get<0>(parsedEqn);
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

int Runfile_c::RunCompare(std::vector<std::complex<mpfr::mpreal>> run1, std::vector<std::complex<mpfr::mpreal>> run2){
	if((run1[0] != run2[0]) || (run1[1] != run2[1]) || (run1[2] != run2[2])) return 0;	// runs are different
	for(unsigned int i = 4; i <= run1.size(); ++i){
		for(unsigned int j = 4; j <= run2.size(); ++j){
			if(run1[i-1] != run2[j-1]) return -2;	// runs differ only by hp
		}
	}
	return -1;							// runs are identical
}
} // namespace virasoro
