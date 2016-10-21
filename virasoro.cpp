#include "virasoro.h"

typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char** argv){
	mpf_set_default_prec(precision);
	if(argc < 6 || !(std::atof(argv[1]) > 0.0) || !(std::atof(argv[2]) > 0) || !(std::atof(argv[3]) > 0) || !(std::atof(argv[4]) >= 0.0) || std::atoi(argv[5]) < 1){
		printf("Error: Central charge c, two dimensions hl and hh, exchange dimension hp, and maximum order N in q^n required as arguments.\n");
		return EXIT_FAILURE;
	}
	auto timeStart = Clock::now();
	if(std::atof(argv[1]) > 1.0 && std::atof(argv[1]) < 25.0){
		printf("Error: Currently not able to accommodate central charge between 1 and 25 due to intermediate complex numbers. Please try again later.\n");
		return EXIT_FAILURE;
	}
	mpf_t c, hl, hh, hp, temp1, temp2;
	mpf_inits(temp1, temp2, NULL);
	mpf_init_set_str(c, argv[1], 10);
	mpf_init_set_str(hl, argv[2], 10);
	mpf_init_set_str(hh, argv[3], 10);
	mpf_init_set_str(hp, argv[4], 10);		
	if(std::atoi(argv[5]) > 65535){
		std::cout << "Error: this version of the program can not accommodate orders higher than 65,535. To remove this restriction, delete all instances of \"short\" and the statement which produces this message." << std::endl;
		return EXIT_FAILURE;
	}
	unsigned short int maxOrder = std::atoi(argv[5]);
	maxOrder -= (maxOrder % 2);
	powOverflow = new mpf_t[maxOrder/256+1];
	mpf_init_set_ui(powOverflow[0], 1);
	for(int i = 1; i <= maxOrder/256; ++i){
		mpf_init_set_ui(powOverflow[i], 16);
		mpf_pow_ui(powOverflow[i], powOverflow[i], 256*i);
	}
	
	int mnLocation[maxOrder]; 	/* "pos" (location+1) in mn vector at which i+1 = m*n/2 starts */
	int mnMultiplicity[maxOrder] = {0};	/* number of mn combinations giving i+1 = m*n/2 */
	int numberOfMN = EnumerateMN(mnLocation, mnMultiplicity, maxOrder);
	unsigned short int mTable[numberOfMN] = {0};
	unsigned short int nTable[numberOfMN] = {0};
	int mnLookup[maxOrder*maxOrder];
	FillMNTable(mnLookup, mTable, nTable, numberOfMN, mnLocation, mnMultiplicity, maxOrder);
	
	// construct b^2 and 1/b^2 from c and lambda_h and lambda_l from h_h and h_l
	mpf_t bsq, invBsq, hpmn[numberOfMN], llsq, lhsq;
	mpf_inits(bsq, invBsq, llsq, lhsq, NULL);
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
	
	FillHpmn(hpmn, mTable, nTable, numberOfMN, bsq, invBsq);
	
	mpf_t prodLkl[numberOfMN];
	auto time1 = Clock::now();
	FillProdLkl(prodLkl, bsq, invBsq, mTable, nTable, numberOfMN);
	auto time2 = Clock::now();
	std::cout << "ProdLkl filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// combine prodLKL into Rmn
	mpf_t Rmn[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos) mpf_init(Rmn[pos-1]);
	time1 = Clock::now();
	FillRmn(Rmn, prodLkl, bsq, invBsq, llsq, lhsq, numberOfMN, mnLookup, maxOrder);
	time2 = Clock::now();
	std::cout << "Rmn filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// combine Rmn and hpmn into computation of H
	mpf_t* Hmn[maxOrder/2];
	Hmn[0] = new mpf_t[numberOfMN];
	for(int pos = 1; pos <= numberOfMN; ++pos) mpf_init_set_ui(Hmn[0][pos-1], 1);
	int mnAtThisOrder;
	for(int order = 2; order <= maxOrder-2; order+=2){
		mnAtThisOrder = 0;
		for(int mn = 2; mn <= maxOrder-order; mn+=2){
			mnAtThisOrder += mnMultiplicity[mn/2-1];
		}
//		std::cout << "Creating " << mnAtThisOrder << " Hmn at order " << order << "." << std::endl;
		Hmn[order/2] = new mpf_t[mnAtThisOrder];
		for(int pos = 1; pos <= mnAtThisOrder; ++pos) mpf_init(Hmn[order/2][pos-1]);
	}
	time1 = Clock::now();
	FillHmn(Hmn, Rmn, hpmn, mnLocation, mnMultiplicity, maxOrder);
	time2 = Clock::now();
	std::cout << "Hmn filled. This took " << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << " us." << std::endl;
	
	// corral H by q order and display coefficients
	mpf_t H[maxOrder/2+1];
	mpf_init_set_ui(H[0], 1);
	for(int i = 1; i <= maxOrder/2; ++i) mpf_init(H[i]);
	FillH(H, Hmn, Rmn, hpmn, hp, mnLocation, mnMultiplicity, maxOrder);
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
	DisplayH(H, c, hl, hh, hp, maxOrder, elapsed, unit);
	
	return EXIT_SUCCESS;
}

int EnumerateMN (int* mnLocation, int* mnMultiplicity, const unsigned short int maxOrder){
	int numberOfMN = 0;
	for(int m=1; m <= maxOrder; m+=2){
		for(int n=2; m*n <= maxOrder; n+=2){		// odd m, even n
			++mnMultiplicity[m*n/2-1];
			++numberOfMN;
		}
		for(int n=1; (m+1)*n <= maxOrder; ++n){		// even m
			++mnMultiplicity[(m+1)/2*n-1];
			++numberOfMN;
		}
	}
	mnLocation[0] = 1;
	for(int i = 2; i <= maxOrder; ++i){
		mnLocation[i-1] = mnLocation[i-2] + mnMultiplicity[i-2];
	}
	return numberOfMN;
}

// currently not able to deal with imaginary bsq from 1 < c < 25.
/*__float128 CToB (const __float128 c){
	__float128 bsq = 0;
	if(c >= 1.0 && c <= 25.0){
		// bad complex value
	} else {
		bsq = c - 13.0q - sqrtq(c*c - 26.0q*c + 25.0q);
		bsq /= 12.0q;
	}
	return bsq;
}*/

void FillMNTable (int *mnLookup, unsigned short int *mTable, unsigned short int *nTable, const int numberOfMN, const int *mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	int pos;
	for(int m = 1; m <= maxOrder; m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){		// odd m, even n
			for(pos = mnLocation[m*n/2-1]; pos <= mnLocation[m*n/2-1]+mnMultiplicity[m*n/2-1]; ++pos){
				if(mTable[pos-1]==0){
					mTable[pos-1] = m;
					nTable[pos-1] = n;
					mnLookup[(m-1)*maxOrder + n-1] = pos;	
					break;
				}
			}
		}
		for(int n = 1; (m+1)*n <= maxOrder; ++n){	// even m
			for(pos = mnLocation[(m+1)*n/2-1]; pos <= mnLocation[(m+1)*n/2-1]+mnMultiplicity[(m+1)*n/2-1]; ++pos){
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

void FillHpmn (mpf_t* hpmn, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN, const mpf_t bsq, const mpf_t invBsq){
	mpf_t temp;
	mpf_init(temp);
	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_init(hpmn[pos-1]);
		mpf_set_ui(hpmn[pos-1], 1);
		mpf_sub_ui(hpmn[pos-1], hpmn[pos-1], nTable[pos-1]*nTable[pos-1]);
		mpf_mul(hpmn[pos-1], hpmn[pos-1], bsq);
		mpf_div_ui(hpmn[pos-1], hpmn[pos-1], 4);
		mpf_set_ui(temp, 1);
		mpf_sub_ui(temp, temp, mTable[pos-1]*mTable[pos-1]);
		mpf_mul(temp, temp, invBsq);
		mpf_div_ui(temp, temp, 4);
		mpf_add(hpmn[pos-1], hpmn[pos-1], temp);
		mpf_set_ui(temp, 1);
		mpf_div_ui(temp, temp, 2);
		mpf_add(hpmn[pos-1], hpmn[pos-1], temp);
		mpf_set_ui(temp, mTable[pos-1]*nTable[pos-1]);
		mpf_div_ui(temp, temp, 2);
		mpf_sub(hpmn[pos-1], hpmn[pos-1], temp);
	}
	mpf_clear(temp);
	return;
}

void FillProdLkl(mpf_t *prodLkl, const mpf_t bsq, const mpf_t invBsq, const unsigned short int *mTable, const unsigned short int *nTable, const int numberOfMN){
	mpf_t temp1, temp2, prod;
	mpf_inits(temp1, temp2, prod, NULL);
	for(int pos = 1; pos <= numberOfMN; ++pos){
		FindProdLkl(bsq, invBsq, mTable[pos-1], nTable[pos-1], temp1, temp2, prod);
		mpf_init_set(prodLkl[pos-1], prod);
	}
	mpf_clears(temp1, temp2, prod, NULL);
}

void FindProdLkl(const mpf_t bsq, const mpf_t invBsq, const unsigned int m, const unsigned int n, mpf_t& temp1, mpf_t& temp2, mpf_t& prod){
//	bool verbose = false;
//	if(m == 2 && n == 3) verbose = true;
	mpf_set_ui(prod,1);
	for(unsigned int l = 1; l <= n-1; ++l){ 		// pairs of Lml*Lm(-l)
		mpf_set_ui(temp1, m*m);
		mpf_mul(temp1, temp1, invBsq);
		mpf_set_ui(temp2, l*l);
		mpf_mul(temp2, temp2, bsq);
		mpf_sub(temp1, temp1, temp2);
		mpf_mul(prod, prod, temp1);
//		if(verbose) std::cout << "Lml*Lm(-l) term." << std::endl;		
	}
	for(unsigned int k = 1; k <= m-1; ++k){ 		// pairs of Lkn*L(-k)n
		mpf_set_ui(temp1, k*k);
		mpf_mul(temp1, temp1, invBsq);
		mpf_set_ui(temp2, n*n);
		mpf_mul(temp2, temp2, bsq);
		mpf_sub(temp1, temp2, temp1);
		mpf_mul(prod, prod, temp1);
//		if(verbose) std::cout << "Lkn*L(-k)n term." << std::endl;		
	}
	for(unsigned int k = 1; k < m; ++k){ 		// pairs of Lk0*L(-k)0
		mpf_set_si(temp1, k*k);
		mpf_mul(temp1, temp1, invBsq);
		mpf_neg(temp1, temp1);
		mpf_mul(prod, prod, temp1);		
//		if(verbose) std::cout << "Lk0*L(-k)0 term." << std::endl;			
	}
	for(unsigned int l = 1; l < n; ++l){ 		// pairs of L0l*L0(-l)
		mpf_set_ui(temp1, l*l);
		mpf_mul(temp1, temp1, bsq);
		mpf_neg(temp1, temp1);
		mpf_mul(prod, prod, temp1);
//		if(verbose) std::cout << "L0l*L0(-l) term." << std::endl;		
	}
	mpf_mul_ui(prod, prod, m*n);
	mpf_neg(prod, prod);				// loose Lm0 and L0n
	for(unsigned int k = 1; k <= m-1; ++k){		
		for(unsigned int l = 1; l <= n-1; ++l){
			mpf_set_ui(temp1, k*k);
			mpf_mul(temp1, temp1, invBsq);
			mpf_set_ui(temp2, l*l);
			mpf_mul(temp2, temp2, bsq);
			mpf_add(temp1, temp1, temp2);
			mpf_add_ui(temp2, temp1, 2*k*l);		// paired L(-k)l*Lk(-l)
			mpf_sub_ui(temp1, temp1, 2*k*l);		// paired Lkl*L(-k)(-l) 
			mpf_mul(prod, prod, temp1);
			mpf_mul(prod, prod, temp2);
//			if(verbose) std::cout << "L" << k << l << "*L(-" << k << ")(-" << l << ") term." << std::endl;
//			if(verbose) std::cout << "L(-" << k << ")" << l << "*L" << k << "(-" << l << ") term." << std::endl;						
		}
	}
	mpf_ui_div(prod, 1, prod);
	return;
}

void FillRmn(mpf_t* Rmn, const mpf_t* prodLkl, const mpf_t bsq, const mpf_t invBsq, const mpf_t llsq, const mpf_t lhsq, const int numberOfMN, const int* mnLookup, const unsigned short int maxOrder){
	mpf_t temp1, temp2, Lsq;
	mpf_inits(temp1, temp2, Lsq, NULL);
	
	mpf_div_ui(temp1, bsq, 16);
	mpf_set(temp2, temp1);
	mpf_add(temp1, temp1, llsq);
	mpf_add(temp2, temp2, lhsq);
	mpf_mul(Rmn[0], temp1, temp2);
	mpf_mul(Rmn[0], Rmn[0], bsq);
	mpf_mul(Rmn[0], Rmn[0], bsq);
	mpf_div_ui(Rmn[0], Rmn[0], 2);
	mpf_neg(Rmn[0], Rmn[0]);																	// R12

	mpf_div_ui(temp1, invBsq, 16);
	mpf_set(temp2, temp1);
	mpf_add(temp1, temp1, llsq);
	mpf_add(temp2, temp2, lhsq);
	mpf_mul(Rmn[1], temp1, temp2);
	mpf_mul(Rmn[1], Rmn[1], invBsq);
	mpf_mul(Rmn[1], Rmn[1], invBsq);
	mpf_div_ui(Rmn[1], Rmn[1], 2);
	mpf_neg(Rmn[1], Rmn[1]);																	// R21


//	Rmn[3] = (0.0625q*(2.0q + invBsq + bsq) + llsq)*(0.0625q*(2.0q + invBsq + bsq) + lhsq)*(2.0q + invBsq + bsq);	// R22
//	Rmn[3] *= (0.0625q*(-2.0q + invBsq + bsq) + llsq)*(0.0625q*(-2.0q + invBsq + bsq) + lhsq)*(-2.0q + invBsq + bsq);
	mpf_set_ui(temp1, 2);
	mpf_add(temp1, temp1, invBsq);
	mpf_add(temp1, temp1, bsq);
	mpf_div_ui(temp1, temp1, 16);
	mpf_add(temp1, temp1, llsq);
	mpf_set_ui(temp2, 2);
	mpf_add(temp2, temp2, invBsq);
	mpf_add(temp2, temp2, bsq);
	mpf_div_ui(temp2, temp2, 16);
	mpf_add(temp2, temp2, lhsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_set_ui(temp2, 2);
	mpf_add(temp2, temp2, invBsq);
	mpf_add(temp2, temp2, bsq);
	mpf_mul(Rmn[3], temp1, temp2);
	mpf_mul(Rmn[3], Rmn[3], temp2);
	mpf_set_ui(temp1, 2);
	mpf_neg(temp1, temp1);
	mpf_add(temp1, temp1, invBsq);
	mpf_add(temp1, temp1, bsq);
	mpf_div_ui(temp1, temp1, 16);
	mpf_add(temp1, temp1, llsq);
	mpf_set_ui(temp2, 2);
	mpf_neg(temp2, temp2);
	mpf_add(temp2, temp2, invBsq);
	mpf_add(temp2, temp2, bsq);
	mpf_div_ui(temp2, temp2, 16);
	mpf_add(temp2, temp2, lhsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_set_ui(temp2, 2);
	mpf_neg(temp2, temp2);
	mpf_add(temp2, temp2, invBsq);
	mpf_add(temp2, temp2, bsq);
	mpf_mul(temp1, temp1, temp2);
	mpf_mul(Rmn[3], Rmn[3], temp1);
	mpf_mul(Rmn[3], Rmn[3], temp2);
	mpf_div_ui(Rmn[3], Rmn[3], 2);
	mpf_neg(Rmn[3], Rmn[3]);											// R22
	
	for(int m = 1; m <= 2; ++m){
		for(int n = 3+(m%2); m*n <= maxOrder; n+=(1+m%2)){
			mpf_set(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-3]-1]);
			for(int p = -m+1; p <= m-1; p+=2){
				mpf_set_ui(Lsq, 2*n-1);
				mpf_sub_ui(Lsq, Lsq, n*n);
				mpf_mul(Lsq, Lsq, bsq);
				mpf_set_si(temp1, 2*(n-1)*p);
				mpf_add(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, invBsq, p*p);
				mpf_sub(Lsq, Lsq, temp1);
//				mpf_set_d(Lsq, (-1.0 + 2.0*n - n*n)*_bsq + 2.0*p*(1-n) - p*p*_invBsq);
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, llsq);
				mpf_sub(temp2, temp2, lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], temp1);
			}
		}
	}

	for(int m = 3; m <= maxOrder-(maxOrder%2); m+=2){
		for(int n = 2; m*n <= maxOrder; n+=2){
			mpf_set(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-3)*maxOrder + n-1]-1]);
			for(int q = -n+1; q <= n-1; q+=2){
				mpf_set_ui(Lsq, 2*m-1);
				mpf_sub_ui(Lsq, Lsq, m*m);
				mpf_mul(Lsq, Lsq, invBsq);
				mpf_set_si(temp1, 2*(m-1)*q);
				mpf_add(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, bsq, q*q);
				mpf_sub(Lsq, Lsq, temp1);			
//				mpf_set_d(Lsq, (-1.0 + 2.0*m - m*m)*_invBsq + 2.0*q*(1-m) - q*q*_bsq);
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, llsq);
				mpf_sub(temp2, temp2, lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], Rmn[mnLookup[(m-1)*maxOrder + n-1]-1], temp1);
			}
		}
		for(int n = 1; m*n+n <= maxOrder; ++n){
			mpf_set(Rmn[mnLookup[m*maxOrder + n-1]-1], Rmn[mnLookup[(m-2)*maxOrder + n-1]-1]);
			for(int q = -n+1; q <= n-1; q+=2){
				mpf_mul_ui(Lsq, invBsq, m*m);
				mpf_neg(Lsq, Lsq);
				mpf_set_si(temp1, 2*q*m);
				mpf_sub(Lsq, Lsq, temp1);
				mpf_mul_ui(temp1, bsq, q*q);
				mpf_sub(Lsq, Lsq, temp1);
//				mpf_set_d(Lsq, -m*m*_invBsq - 2.0*q*m - q*q*_bsq);
				mpf_div_ui(temp1, Lsq, 16);
				mpf_set(temp2, temp1);
				mpf_sub(temp1, temp1, llsq);
				mpf_sub(temp2, temp2, lhsq);
				mpf_mul(temp1, temp1, temp2);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(temp1, temp1, Lsq);
				mpf_mul(Rmn[mnLookup[m*maxOrder + n-1]-1], Rmn[mnLookup[m*maxOrder + n-1]-1], temp1);
			}
		}
	}

	for(int pos = 1; pos <= numberOfMN; ++pos){
		mpf_mul(Rmn[pos-1], Rmn[pos-1], prodLkl[pos-1]);
		std::cout << "Rmn[" << pos-1 << "] = ";
		mpf_out_str(NULL, 10, 10, Rmn[pos-1]);
		std::cout << "." << std::endl;		
	}
	
	mpf_clears(temp1, temp2, Lsq);
	return;
}

void FillHmn(mpf_t** Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	std::thread thread[maxThreads];
	mpf_t temp1[maxThreads], temp2[maxThreads], hpTemp[maxThreads];
	for(int i=1; i <= maxThreads; ++i) mpf_inits(temp1[i-1], temp2[i-1], hpTemp[i-1], NULL);
	int numThreads;
	for(int order = 2; order < maxOrder; order+=2){
		numThreads = std::min(maxThreads,(maxOrder-order)/2);
//		std::cout << "Launching thread to fill mn = " << 2 << " through " << (maxOrder-order)/numThreads + (maxOrder-order)%numThreads << " of " << mnLocation[maxOrder/2-1] + mnMultiplicity[maxOrder/2-1] << "." << std::endl;
		thread[0] = std::thread(ThreadFillHmn, Hmn, Rmn, hpmn, mnLocation, mnMultiplicity, 2, (maxOrder-order)/numThreads + (maxOrder-order)%numThreads, order, std::ref(temp1[0]), std::ref(temp2[0]), std::ref(hpTemp[0]));
		for(int i=2; i<=numThreads; ++i){
//			std::cout << "Launching thread to fill mn = " << (i-1)*(maxOrder-order)/numThreads + (maxOrder-order)%numThreads + 1 << " through " << i*(maxOrder-order)/numThreads + (maxOrder-order)%numThreads << " of " << mnLocation[maxOrder/2-1] + mnMultiplicity[maxOrder/2-1] << "." << std::endl;		
			thread[i-1] = std::thread(ThreadFillHmn, Hmn, Rmn, hpmn, mnLocation, mnMultiplicity, (i-1)*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads + 1, i*((maxOrder-order)/numThreads) + (maxOrder-order)%numThreads, order, std::ref(temp1[i-1]), std::ref(temp2[i-1]), std::ref(hpTemp[i-1]));
		}
		for(int i=1; i<= numThreads; ++i){
			thread[i-1].join();
		}
//		std::cout << "--- ** Filled order " << order << ". ** ---" << std::endl;
	}
//	std::cout << "Last Hmn filled." << std::endl;
	for(int i=1; i <= maxThreads; ++i) mpf_clears(temp1[i-1], temp2[i-1], hpTemp[i-1], NULL);
	return;
}

void ThreadFillHmn(mpf_t** Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const int* mnLocation, const int* mnMultiplicity, const int startingMN, const int endingMN, const int order, mpf_t& temp1, mpf_t& temp2, mpf_t& hpTemp){
	for(int mn = startingMN + startingMN%2; mn <= endingMN; mn+=2){
//		std::cout << "Filling mn = " << mn << " at order " << order << ". I think there are " << mnMultiplicity[mn/2-1] << " Hmn here, with positions " << mnLocation[mn/2-1]-1 << " through " << mnLocation[mn/2-1]+mnMultiplicity[mn/2-1]-2 << "." << std::endl;
		for(int pos = mnLocation[mn/2-1]; pos <= mnLocation[mn/2-1] + mnMultiplicity[mn/2-1] - 1; ++pos){
//			std::cout << "Attempting to fill Hmn[" << pos-1 << "] on order " << order << "." << std::endl;
			mpf_add_ui(hpTemp, hpmn[pos-1], mn);
			for(int power = 2; power <= order; power+=2){
				for(int scanPos = mnLocation[power/2-1]; scanPos <= mnLocation[power/2-1] + mnMultiplicity[power/2-1] - 1; ++scanPos){
//					std::cout << "Finding contribution to Hmn[" << order << "," << pos-1 << "] from Hmn[" << (order-power) << "," << scanPos-1 << "]." << std::endl;
					mpf_sub(temp1, hpTemp, hpmn[scanPos-1]);
					mpf_div(temp1, Rmn[scanPos-1], temp1);
					mpf_mul(temp1, temp1, Hmn[(order-power)/2][scanPos-1]);
					mpf_mul(temp1, temp1, powOverflow[power/256]);
					mpf_set_ui(temp2, 16);
					mpf_pow_ui(temp2, temp2, power%256);
					mpf_mul(temp1, temp1, temp2);
					mpf_add(Hmn[order/2][pos-1], Hmn[order/2][pos-1], temp1);
				}
			}
//			std::cout << "Filled Hmn[" << pos-1 << "] on order " << order << "." << std::endl;
		}
	}
	return;
}

void FillH(mpf_t* H, const mpf_t* const* Hmn, const mpf_t* Rmn, const mpf_t* hpmn, const mpf_t hp, const int* mnLocation, const int* mnMultiplicity, const unsigned short int maxOrder){
	mpf_t temp1, temp2;
	mpf_inits(temp1, temp2, NULL);
	for(int order = 2; order <= maxOrder; order+=2){
		for(int power = 2; power <= order; power+=2){
			for(int scanPos = mnLocation[power/2-1]; scanPos <= mnLocation[power/2-1] + mnMultiplicity[power/2-1] - 1; ++scanPos){
				mpf_sub(temp1, hp, hpmn[scanPos-1]);
				mpf_div(temp1, Rmn[scanPos-1], temp1);
				mpf_mul(temp1, temp1, Hmn[(order-power)/2][scanPos-1]);
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

inline std::string to_string(const mpf_t N, const int digits){
	char* buffer = new char;
	mp_exp_t* dotPos = new long;
	buffer = mpf_get_str(NULL, dotPos, 10, digits, N);
	std::string output = std::string(buffer);
	if(output.empty()) output.append("0");
	while((long)*dotPos > (int)output.size()) output.append("0");
	if((long)*dotPos > 0 && (long)*dotPos < (int)output.size()) output.insert((long)*dotPos, ".");
	delete buffer;
	delete dotPos;
	return output;
}

void DisplayH(const mpf_t* H, const mpf_t c, const mpf_t hl, const mpf_t hh, const mpf_t hp, const unsigned short int maxOrder, const int time, const std::string unit){
	std::ofstream outputFile;
	std::string filename = "virasoro_" + to_string(c, 1) + "_" + to_string(hl, 1) + "_" + to_string(hh, 1) + "_" + to_string(hp, 1) + "_" + std::to_string(maxOrder) + ".txt";
	outputFile.open (filename);
	std::cout << "Given the parameters" << std::endl;
	outputFile << "Given the parameters" << std::endl;
	std::cout << "c = " << to_string(c, 10) << ", h_L = " << to_string(hl, 10) << ", h_H = " << to_string(hh, 10) << ", h_p = " << to_string(hp, 10) << std::endl;
	outputFile << "c = " << to_string(c, 0) << ", h_L = " << to_string(hl, 0) << ", h_H = " << to_string(hh, 0) << ", h_p = " << to_string(hp, 0) << std::endl;	
	std::cout << "the Virasoro block coefficients are as follows:" << std::endl;
	outputFile << "the Virasoro block coefficients are as follows:" << std::endl;
	for(int orderBy2 = 0; 2*orderBy2 <= maxOrder; ++orderBy2){
		std::cout << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 10) << std::endl;
		outputFile << "q^" << 2*orderBy2 << ": " << to_string(H[orderBy2], 0) << std::endl;
	}
	outputFile << "{1";
	for(int orderBy2 = 1; 2*orderBy2 <= maxOrder; orderBy2++){
		outputFile << "," << to_string(H[orderBy2], 0);
	}
	outputFile << "}" << std::endl;
	outputFile.close();
	std::cout << "This computation took " << time << unit << "." << std::endl;
}
