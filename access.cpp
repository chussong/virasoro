#include "access.h"

namespace virasoro {
int Access::maxOrder = 0;
int Access::numberOfMN = 0;
std::vector<int> Access::mTable; // mTable[loc] = m at loc
std::vector<int> Access::nTable; // nTable[loc] = n at loc
std::vector<int> Access::mnLookup; // mnLookup[(m-1)*maxOrder + n-1] = pos of (m,n)
std::vector<int> Access::mnLocation; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
std::vector<int> Access::mnMultiplicity;	/* number of mn combinations giving i+1 = m*n */

void Access::Populate(const int maxOrder){
	if(Access::maxOrder >= maxOrder) return;
	mnLocation.resize(maxOrder);
	mnMultiplicity.resize(maxOrder);
	for(int m=1; m <= maxOrder; m+=2){
		for(int n=std::max(2,Access::maxOrder/m+2); m*n <= maxOrder; n+=2-(n%2)){			// odd m, even n
			++mnMultiplicity[m*n-1];
			++numberOfMN;
		}
		for(int n=std::max(1,Access::maxOrder/(m+1)+1); (m+1)*n <= maxOrder; ++n){	// even m
			++mnMultiplicity[(m+1)*n-1];
			++numberOfMN;
		}
	}
	mnLocation[1] = 1;
	for(int i = std::max(4,Access::maxOrder+2); i <= maxOrder; i+=2){
		mnLocation[i-1] = mnLocation[i-3] + mnMultiplicity[i-3];
	}
	mTable.resize(numberOfMN);
	nTable.resize(numberOfMN);
	mnLookup.resize(maxOrder*maxOrder);
	int pos;
	for(int m = 1; m <= maxOrder; m+=2){
		for(int n = std::max(2,Access::maxOrder/m+2); m*n <= maxOrder; n+=2-(n%2)){		// odd m, even n
			for(pos = mnLocation[m*n-1]; pos <= mnLocation[m*n-1]+mnMultiplicity[m*n-1]; ++pos){
				if(mTable[pos-1]==0){
					mTable[pos-1] = m;
					nTable[pos-1] = n;
					mnLookup[(m-1)*maxOrder + n-1] = pos;	
					break;
				}
			}
		}
		for(int n = std::max(1,Access::maxOrder/(m+1)+1); (m+1)*n <= maxOrder; ++n){	// even m, all n
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
	Access::maxOrder = maxOrder;
}

int Access::MNAtLevel(const int level){
	if(level % 2 == 1 || level < 0 || level > maxOrder) return 0;
	int ret = 0;
	for(int mn = 2; mn < maxOrder; mn+=2){
		ret += MultOfMN(mn);
	}
	return ret;
}

// digits <= 0 gives exact number
std::string to_string(const mpfr::mpreal& N, int digits){
	if(digits <= 0) return N.toString(-1, 10, MPFR_RNDN);
	std::string output = N.toString(digits, 10, MPFR_RNDN);
	std::size_t ePos = output.find('e');
	if(ePos < std::string::npos) output = output.replace(ePos, 1, "*10^");
	return output;
}

// digits = 0 gives exact number in Mathematica syntax, digits < 0 gives exact
// number in MPC syntax
std::string to_string(const std::complex<mpfr::mpreal>& N, int digits, int base){
	char* cstr = mpc_get_str(base, std::max(digits,0), N.mpc_srcptr(), MPC_RNDNN);
	std::string output(cstr);
	mpc_free_str(cstr);
	if(digits < 0){
		return output;
	}
	size_t splitLoc = output.find(" ");
	std::string halves[2];
	halves[0] = output.substr(1, splitLoc-1);
	halves[1] = output.substr(splitLoc+1, output.size()-splitLoc-2);
	mpfr::mpreal mpfHalf;
	for(int i = 1; i <= 2; ++i){
		if(halves[i-1] == "+0" || halves[i-1] == "-0"){
			halves[i-1].clear();
		} else {
			if(halves[i-1][0] == '-'){
				mpfHalf = halves[i-1].substr(1);
			} else {
				mpfHalf = halves[i-1];
			}
			if(mpfHalf < STATICTOLERANCE) halves[i-1].clear(); // change to regular tolerance
		}
	}
	size_t eLoc, eEnd;
	for(int i = 1; i <= 2; ++i){
		if(halves[i-1].empty()) continue;
		if((eLoc=halves[i-1].find("e")) < std::string::npos){
			eEnd = halves[i-1].find_first_of(" )", eLoc+3);
			int exp = std::stoi(halves[i-1].substr(eLoc+1, eEnd-eLoc-1));
			if(digits == 0 || std::abs(exp) < digits){
				halves[i-1].erase(eLoc, eEnd-eLoc);
				if(exp > 0){
					if(halves[i-1][0] == '-'){
						halves[i-1].erase(2, 1);
						if(static_cast<unsigned long>(exp)+2 >= halves[i-1].size()) halves[i-1].append(exp+3-halves[i-1].size(), '0');
						halves[i-1].insert(exp+2, ".");
					} else {
						halves[i-1].erase(1, 1);
						if(static_cast<unsigned long>(exp)+1 >= halves[i-1].size()) halves[i-1].append(exp+2-halves[i-1].size(), '0');
						halves[i-1].insert(exp+1, ".");
					}
				} else {
					if(halves[i-1][0] == '-'){
						halves[i-1].erase(2, 1);
						halves[i-1].insert(1, "0.");
						halves[i-1].insert(3, -1-exp, '0');
					} else {
						halves[i-1].erase(1, 1);
						halves[i-1].insert(0, "0.");
						halves[i-1].insert(2, -1-exp, '0');
					}
				}
			} else {
				if(halves[i-1][eLoc+1] == '+'){
					halves[i-1].replace(eLoc, 2, "*10^");
				} else {
					halves[i-1].replace(eLoc, 1, "*10^");
				}
			}
		}
		eEnd = halves[i-1].find_last_not_of("0");
		if(halves[i-1][eEnd] == '.') halves[i-1].erase(eEnd);
	}
	if(halves[0].empty()){
		if(halves[1].empty()){
			return "0";
		} else {
			halves[1].append("*I");
			return halves[1];
		}
	} else {
		if(halves[1].empty()){
			return halves[0];
		} else {
			if(halves[1][0] == '-') return halves[0] + "-" + halves[1].substr(1) + "*I";
			return halves[0] + "+" + halves[1] + "*I";
		}
	}
}

} // namespace virasoro
