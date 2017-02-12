#include "access.h"

//namespace virasoro {
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
//} // namespace virasoro
