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
	entries = ParseLines(lines);
	int errorCode = Expand();
	if(errorCode == -2){
		perror("Error: a curly brace '{' was found indicating a batch run, but it was not followed by a valid macro.\n");
		return -2;
	}
	std::tuple<std::vector<std::vector<std::complex<mpfr::mpreal>>>,std::vector<int>> numericization = NumericizeRuns(entries);
	runs = std::get<0>(numericization);
	maxOrders = std::get<1>(numericization);
	if(runs.empty()){
		perror("Error: zero valid runs detected in given runfile.\n");
		return 0;
	}
	CrunchDuplicateRuns(runs, maxOrders);
	lines.clear();
	entries.clear();
	return runs.size();
}

std::vector<std::vector<std::string>> Runfile_c::ParseLines(const std::vector<std::string>& lines){
	std::vector<std::vector<std::string>> entries;
	size_t gapPos = -1;
	size_t entryStart;
	std::vector<std::string> entryVec;
	for(unsigned int i = 0; i < lines.size(); ++i){
		do{
			entryStart = lines[i].find_first_not_of(' ', gapPos+1);
			if(entryStart == std::string::npos) break;
			gapPos = lines[i].find_first_of(' ', entryStart+1);
			if(lines[i].find_first_of('(', entryStart) < gapPos){
				gapPos = lines[i].find_first_of(' ', gapPos+1);
			}
			entryVec.push_back(lines[i].substr(entryStart, gapPos-entryStart));
		}while(gapPos != std::string::npos);
		entries.push_back(entryVec);
		entryVec.clear();
	}
	return entries;
}

int Runfile_c::Expand(){
	for(int param = 0; param < 4; ++param){
		ExpandRelativeEqns(param);
		ExpandBraces(param);
	}
	return 0; 
}

// param=0 is c/b/b^2; param=1 is hl; param=2 is hh; param=3 is hp
int Runfile_c::ExpandBraces(const int param){
	std::vector<std::vector<std::string>> newEntries;
	std::tuple<std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>, std::complex<mpfr::mpreal>> parsedBraces;
	std::complex<mpfr::mpreal> currentValue;
	std::string insideBraces;
	size_t leftPos, rightPos;
	for(unsigned int i = 0; i < entries.size(); ++i){
		if((leftPos = entries[i][param].find("{")) != std::string::npos){
			if((rightPos = entries[i][param].find("}", leftPos+1)) == std::string::npos){
				return -2;
			} else if(entries[i][param].find_first_not_of("0123456789-+()e. ,;", leftPos+1) == rightPos) {
				insideBraces = entries[i][param].substr(leftPos+1, rightPos-leftPos-1);
				parsedBraces = ParseBraces(insideBraces);
				for(currentValue = std::get<0>(parsedBraces); currentValue.real() <= std::get<1>(parsedBraces).real(); currentValue += std::get<2>(parsedBraces)){
					newEntries.push_back(entries[i]);
					newEntries[newEntries.size()-1][param] = to_string(currentValue, -1);
				}
			} else {
				newEntries.push_back(entries[i]);
			}
		} else {
			newEntries.push_back(entries[i]);
		}
	}
	entries = newEntries;
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

// param=0 is c/b/b^2; param=1 is hl; param=2 is hh; param=3 is hp
void Runfile_c::ExpandRelativeEqns(const int param){
	std::vector<std::string> newEntries;
	std::string firstHalf, secondHalf;
	size_t leftPos, rightPos, splitPos;
	std::complex<mpfr::mpreal> value;
	bool madeChange = true;
	for(unsigned int i = 0; i < entries.size(); ++i){
		do{
			madeChange = false;
			if((splitPos = entries[i][param].find_first_of("ch")) != std::string::npos){
				leftPos = entries[i][param].find_last_of(" ,;{", splitPos-1) + 1;
				rightPos = entries[i][param].find_first_of(" ,;}", splitPos+1) - 1;
				if(rightPos >= entries[i][param].length()) rightPos = entries[i][param].length()-1;
				firstHalf = entries[i][param].substr(0, leftPos);
				secondHalf = entries[i][param].substr(rightPos+1);
				value = RelativeMPF(i, entries[i][param].substr(leftPos, rightPos-leftPos+1));
				madeChange = true;
				entries[i][param] = firstHalf + to_string(value, -1) + secondHalf;
			}
		}while(madeChange);
	}
	return;
}

std::tuple<std::complex<mpfr::mpreal>, int> Runfile_c::ParseRelativeEqn(std::string equation, std::string relTo){
	std::complex<mpfr::mpreal> modifier(0);
	int type = -100;
	std::size_t hit;
	if((hit = equation.find(relTo)) != std::string::npos){
		if(hit == 0){
			if(equation[relTo.length()] == '+'){			// X + M
				type = 0;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '-'){ 	// X - M
				type = 1;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '*'){		// X * M
				type = 3;
				modifier = equation.substr(relTo.size()+1);
			} else if(equation[relTo.length()] == '/'){		// X / M
				type = 4;
				modifier = equation.substr(relTo.size()+1);
			} else {				// assume it's just equality with no arithmetic
				type = 0;
				modifier = 0;
			}
		} else if(hit >= 2){
			if(equation[hit-1] == '+'){						// M + X
				type = 0;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '-'){				// M - X
				type = 2;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '*'){				// M * X
				type = 3;
				modifier = equation.substr(0,hit-1);
			} else if(equation[hit-1] == '/'){				// M / X
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

std::complex<mpfr::mpreal> Runfile_c::RelativeMPF(int lineNum, std::string equation){
	std::complex<mpfr::mpreal> baseMPF;
	std::tuple<std::complex<mpfr::mpreal>, int> parsedEqn;
	if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "c")) < 0){
		if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hl")) < 0){
			if(std::get<1>(parsedEqn = ParseRelativeEqn(equation, "hh")) < 0){
				return std::complex<mpfr::mpreal>("(+NaN +NaN)");
			}
		}
	}
	std::complex<mpfr::mpreal> modifier = std::get<0>(parsedEqn);
	int type = std::get<1>(parsedEqn);
	baseMPF = entries[lineNum][type/10-1];
	switch(type % 10){
		case 0:		return modifier	+ baseMPF;
		case 1:		return baseMPF	- modifier;
		case 2:		return modifier	- baseMPF;
		case 3:		return modifier	* baseMPF;
		case 4:		return baseMPF	/ modifier;
		case 5:		return modifier	/ baseMPF;
	}
	return std::complex<mpfr::mpreal>("(+NaN +NaN)");
}

std::tuple<std::vector<std::vector<std::complex<mpfr::mpreal>>>,std::vector<int>> Runfile_c::NumericizeRuns(const std::vector<std::vector<std::string>>& entries){
	std::vector<std::vector<std::complex<mpfr::mpreal>>> runs;
	std::vector<std::complex<mpfr::mpreal>> currentRun;
	std::vector<int> maxOrders;
	for(unsigned int i = 0; i < entries.size(); ++i){
		for(int j = 0; j < 4; ++j){
			currentRun.emplace_back(entries[i][j]);
		}
		runs.push_back(currentRun);
		maxOrders.push_back(std::stoi(entries[i][4]));
		currentRun.clear();
	}
	return std::make_tuple(runs, maxOrders);
}

int Runfile_c::RunCompare(const std::vector<std::complex<mpfr::mpreal>>& run1, const std::vector<std::complex<mpfr::mpreal>>& run2){
	if((run1[0] != run2[0]) || (run1[1] != run2[1]) || (run1[2] != run2[2])) return 0;	// runs are different
	for(unsigned int i = 4; i <= run1.size(); ++i){
		for(unsigned int j = 4; j <= run2.size(); ++j){
			if(run1[i-1] != run2[j-1]) return -2;	// runs differ only by hp
		}
	}
	return -1;							// runs are identical
}

// Check runs for duplicates, compress runs differing only by hp
void Runfile_c::CrunchDuplicateRuns(std::vector<std::vector<std::complex<mpfr::mpreal>>>& runs, std::vector<int>& maxOrders){
	std::vector<bool> crunched;
	crunched.resize(runs.size());
	for(unsigned int i = 0; i < runs.size(); ++i){
		for(unsigned int j = i+1; j < runs.size(); ++j){
			if(crunched[j]) continue;
			switch(RunCompare(runs[i], runs[j])){
				case 0:		crunched[j] = false;	// runs different, do nothing
							break;				
				case -1:	if(maxOrders[j] > maxOrders[i]) maxOrders[i] = maxOrders[j];
							crunched[j] = true;	// runs identical, crunch them together
							break;
				case -2:	for(unsigned int k = 3; k < runs[j].size(); ++k) runs[i].push_back(runs[j][k]);
							if(maxOrders[j] > maxOrders[i]) maxOrders[i] = maxOrders[j];
							crunched[j] = true;	// runs differ by hp, make them fast
							break;
			}
		}
	}
	for(int i = runs.size()-1; i >= 0; --i){
		if(crunched[i]){
			runs.erase(runs.begin()+i);
			maxOrders.erase(maxOrders.begin()+i);
		}
	}
	return;
}
} // namespace virasoro
