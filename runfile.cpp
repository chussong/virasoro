#include "runfile.h"

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
	std::vector<mpf_class> currentRun;
	int currentMO;
	for(unsigned int i = 1; i <= lines.size(); ++i){
		leftPos = 0;
		rightPos = 0;
		for(int j = 1; j <= 4; ++j){
			leftPos = lines[i-1].find_first_of("0123456789.-", rightPos);
			rightPos = lines[i-1].find_first_not_of("0123456789.-", leftPos);
			currentRun.emplace_back(lines[i-1].substr(leftPos, rightPos-leftPos));
/*			if(something){
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
	int changesMade;
	int expansions = 1;
	while(expansions > 0){
		expansions = 0;
		changesMade = 1;
		while(changesMade > 0){
			changesMade = ExpandBraces();
			expansions += changesMade;
		}
		changesMade = 1;
		while(changesMade > 0){
			changesMade = ExpandRelativeEqns();
			expansions += changesMade;
		}
	}
	return 0; 
}

int Runfile_c::ExpandBraces(){
	std::vector<std::string> newLines;
	std::string firstHalf, secondHalf, insideBraces;
	std::size_t leftPos, rightPos;
	std::tuple<mpf_class, mpf_class, mpf_class> parsedBraces;
	bool needRerun;
	do{
		needRerun = false;
		for(unsigned int i = 1; i <= lines.size(); ++i){
			if((leftPos = lines[i-1].find("{")) != std::string::npos){
				firstHalf = lines[i-1].substr(0, leftPos);
				if((rightPos = lines[i-1].find("}", leftPos+1)) == std::string::npos){
					return -2;
				} else if(lines[i-1].find_first_not_of("0123456789-. ,;", leftPos+1) == rightPos) {
					secondHalf = lines[i-1].substr(rightPos+1, lines[i-1].length()-rightPos-1);
					insideBraces = lines[i-1].substr(leftPos+1, rightPos-leftPos-1);
					parsedBraces = ParseBraces(firstHalf, insideBraces);
					if(std::get<0>(parsedBraces) > std::get<1>(parsedBraces) || std::get<2>(parsedBraces) <= 0) return -2;
					for(mpf_class currentValue = std::get<0>(parsedBraces); currentValue <= std::get<1>(parsedBraces); currentValue += std::get<2>(parsedBraces)){
						newLines.push_back(firstHalf + to_string(currentValue, 0) + secondHalf);
					}
					needRerun = true;
				} else {
					newLines.push_back(lines[i-1]);
				}
			} else {
				newLines.push_back(lines[i-1]);
			}
		}
		lines = newLines;
		newLines.clear();
	}while(needRerun);
	return 0;
}

std::tuple<mpf_class, mpf_class, mpf_class> Runfile_c::ParseBraces(std::string firstHalf, std::string insideBraces){
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

int Runfile_c::ExpandRelativeEqns(){
	int changesMade = 0;
	std::vector<std::string> newLines;
	std::string currentLine, firstHalf, secondHalf;
	size_t leftPos, rightPos, splitPos;
	mpf_class value;
	bool madeChange = true;
	while(madeChange){
		madeChange = false;
		for(unsigned int i = 1; i <= lines.size(); ++i){
			if((splitPos = lines[i-1].find_first_of("ch")) != std::string::npos){
				leftPos = lines[i-1].find_last_of(" ,;{", splitPos-1);
				rightPos = lines[i-1].find_first_of(" ,;}", splitPos+1);
				firstHalf = lines[i-1].substr(0, leftPos+1);
				secondHalf = lines[i-1].substr(rightPos, currentLine.length() - rightPos);
				value = RelativeMPF(firstHalf, lines[i-1].substr(leftPos+1, rightPos-leftPos-1));
				madeChange = true;
				++changesMade;
				newLines.push_back(firstHalf + to_string(value, 0) + secondHalf);
			} else {
				newLines.push_back(lines[i-1]);
			}
		}
		lines = newLines;
		newLines.clear();
	}
	return changesMade;
}

std::tuple<mpf_class, int> Runfile_c::ParseRelativeEqn(std::string equation, std::string relTo){
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

mpf_class Runfile_c::RelativeMPF(std::string firstHalf, std::string equation){
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
					std::cout << "Dividing " << to_string(baseMPF, 10) << " by " << to_string(modifier, 0) << " to get " << to_string(baseMPF/modifier, 0) << std::endl;
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

int Runfile_c::RunCompare(std::vector<mpf_class> run1, std::vector<mpf_class> run2){
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
