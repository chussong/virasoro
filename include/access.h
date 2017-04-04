#ifndef ACCESS_H_
#define ACCESS_H_
#include <vector>
#include <string>
#include "mpreal.h"
#include "mpcomplex.h"

#define STATICTOLERANCE 10e-100

namespace virasoro {
class Access {
	static int maxOrder;
	static int numberOfMN;
	static std::vector<int> mTable; // mTable[loc] = m at loc
	static std::vector<int> nTable; // nTable[loc] = n at loc
	static std::vector<int> mnLookup; // mnLookup[(m-1)*maxOrder + n-1] = pos of (m,n)
	static std::vector<int> mnLocation; /* "pos" (location+1) in mn vector at which i+1 = m*n starts */
	static std::vector<int> mnMultiplicity;	/* number of mn combinations giving i+1 = m*n */

	public:
		inline static int mAtLoc	(const int loc)				{ return mTable[loc]; }
		inline static int nAtLoc	(const int loc)				{ return nTable[loc]; }
		inline static int mnAtLoc	(const int loc)				{ return mTable[loc]*nTable[loc]; }
		inline static int PosFromMN	(const int m, const int n)	{ return mnLookup[(m-1)*maxOrder + n-1]; }
		inline static int MultOfMN	(const int mn)				{ return mnMultiplicity[mn-1]; }
		inline static int TotalMN	()							{ return numberOfMN;  }

		static void Populate(const int maxOrder);
		static int MNAtLevel(const int level);
};

std::string to_string(const mpfr::mpreal& N, int digits);
std::string to_string(const std::complex<mpfr::mpreal>& N, int digits, int base = 10);

} // namespace virasoro
#endif
