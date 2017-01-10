#ifndef MPCOMPLEX_H_
#define MPCOMPLEX_H_

#include <string>
#include <mpc.h>
#include "mpreal.h"

namespace mpfr{
class mpcomplex{

	public:
		// these two must be defined in a .cpp file (outside of the actual program)
		static mpc_rnd_t default_rnd_mode;
		static mpfr_prec_t default_prec;
		static mpfr::mpreal tolerance;

		static void set_default_prec(const mpfr_prec_t prec);
		static void set_default_rnd_mode(const mpc_rnd_t rnd_mode);
		static void set_tolerance(const mpfr::mpreal tol);

		mpc_rnd_t rnd_mode;
		mpc_t value;

		mpcomplex();
		mpcomplex(const mpcomplex& u);
		mpcomplex(mpcomplex&& u);
		mpcomplex(long rlPart, long imPart = 0, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		mpcomplex(mpfr::mpreal rlPart, mpfr::mpreal imPart = 0, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		mpcomplex(const std::string& str, int base = 10, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		~mpcomplex();
		
		mpcomplex& operator=(const mpcomplex& v);
		mpcomplex& operator=(mpcomplex&& v);
		mpcomplex& operator=(const mpc_t v);
		mpcomplex& operator=(const int v);
		mpcomplex& operator=(const long v);
		mpcomplex& operator=(const unsigned int v);
		mpcomplex& operator=(const unsigned long v);
		mpcomplex& operator=(const char* s);
		mpcomplex& operator=(const std::string& s);

		mpcomplex& operator+=(const mpcomplex& v);
		mpcomplex& operator+=(const mpc_t v);
		mpcomplex& operator+=(const int v);
		mpcomplex& operator+=(const long v);
		mpcomplex& operator+=(const unsigned int v);
		mpcomplex& operator+=(const unsigned long v);

		mpcomplex& operator-=(const mpcomplex& v);
		mpcomplex& operator-=(const mpc_t v);
		mpcomplex& operator-=(const int v);
		mpcomplex& operator-=(const long v);
		mpcomplex& operator-=(const unsigned int v);
		mpcomplex& operator-=(const unsigned long v);
		const mpcomplex operator-() const;

		mpcomplex& operator*=(const mpcomplex& v);
		mpcomplex& operator*=(const mpc_t v);
		mpcomplex& operator*=(const int v);
		mpcomplex& operator*=(const long v);
		mpcomplex& operator*=(const unsigned int v);
		mpcomplex& operator*=(const unsigned long v);

		mpcomplex& operator/=(const mpcomplex& v);
		mpcomplex& operator/=(const mpc_t v);
		mpcomplex& operator/=(const int v);
		mpcomplex& operator/=(const long v);
		mpcomplex& operator/=(const unsigned int v);
		mpcomplex& operator/=(const unsigned long v);

		// <<= is fast multiplication by 2^u
		mpcomplex& operator<<=(const int u);
		mpcomplex& operator<<=(const long u);
		mpcomplex& operator<<=(const unsigned int u);
		mpcomplex& operator<<=(const unsigned long u);

		std::string to_string(int digits = 0, int base = 10) const;
		mpfr::mpreal realPart();
		mpfr::mpreal imPart();
		bool isReal();

		mpcomplex sqr();
		mpcomplex sqrt();

		const mpcomplex sqrt(const mpcomplex& v);
		const mpfr_t* abs(const mpcomplex& v);
		const std::string to_string(const mpcomplex& v, int digits = 0);

		friend void swap(mpcomplex& x, mpcomplex& y);
};

inline void mpcomplex::set_default_prec(const mpfr_prec_t prec){
	default_prec = prec;
}

inline void mpcomplex::set_default_rnd_mode(const mpc_rnd_t rnd_mode){
	default_rnd_mode = rnd_mode;
}

inline void mpcomplex::set_tolerance(const mpfr::mpreal tol){
	tolerance = tol;
}

inline mpcomplex::mpcomplex(){
	rnd_mode = default_rnd_mode;
	mpc_init2(value, mpcomplex::default_prec);
	mpc_set_ui(value, 0, rnd_mode);
}

inline mpcomplex::mpcomplex(const mpcomplex& u){
	rnd_mode = u.rnd_mode;
	mpc_init2(value, mpcomplex::default_prec);
	mpc_set(value, u.value, rnd_mode);
}

inline mpcomplex::mpcomplex(mpcomplex&& u){
	rnd_mode = u.rnd_mode;
	mpc_init2(value, mpcomplex::default_prec);
	mpc_swap(value, u.value);
}

inline mpcomplex::mpcomplex(long rlPart, long imPart, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_si_si(value, rlPart, imPart, rnd_mode);
}

inline mpcomplex::mpcomplex(mpfr::mpreal rlPart, mpfr::mpreal imPart, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_fr_fr(value, rlPart.mpfr_srcptr(), imPart.mpfr_srcptr(), rnd_mode);
}

inline mpcomplex::mpcomplex(const std::string& str, int base, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_str(value, str.c_str(), base, rnd_mode);
}

inline mpcomplex::~mpcomplex(){
	mpc_clear(value);
}

// this is needed for the templates
namespace mpfc_internal{

	template <typename ArgumentType> struct result_type {};

	template <> struct result_type<mpcomplex>			{typedef mpcomplex type;};
	template <> struct result_type<int>					{typedef mpcomplex type;};
	template <> struct result_type<long>				{typedef mpcomplex type;};
	template <> struct result_type<unsigned int>		{typedef mpcomplex type;};
	template <> struct result_type<unsigned long>		{typedef mpcomplex type;};
}

template <typename Rhs>
inline const typename mpfc_internal::result_type<Rhs>::type
	operator+(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) += rhs; }

template <typename Lhs>
inline const typename mpfc_internal::result_type<Lhs>::type
	operator+(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(rhs) += lhs; }

template <typename Rhs>
inline const typename mpfc_internal::result_type<Rhs>::type
	operator-(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) -= rhs; }

template <typename Lhs>
inline const typename mpfc_internal::result_type<Lhs>::type
	operator-(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(lhs) -= rhs; }

template <typename Rhs>
inline const typename mpfc_internal::result_type<Rhs>::type
	operator*(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) *= rhs; }

template <typename Lhs>
inline const typename mpfc_internal::result_type<Lhs>::type
	operator*(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(rhs) *= lhs; }

template <typename Rhs>
inline const typename mpfc_internal::result_type<Rhs>::type
	operator/(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) /= rhs; }

template <typename Lhs>
inline const typename mpfc_internal::result_type<Lhs>::type
	operator/(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(lhs) /= rhs; }

inline mpcomplex& mpcomplex::operator=(const mpcomplex& v){
	mpc_set(value, v.value, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(mpcomplex&& v){
	mpc_swap(value, v.value);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpc_t v){
	mpc_set(value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const int v){
	mpc_set_si(value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const long v){
	mpc_set_si(value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const unsigned int v){
	mpc_set_ui(value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const unsigned long v){
	mpc_set_ui(value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const char* s){
	mpc_set_str(value, s, 10, rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const std::string& s){
	mpc_set_str(value, s.c_str(), 10, rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpcomplex& v){
	mpc_add(value, value, v.value, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpc_t v){
	mpc_add(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const int v){
	mpc_add_si(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const long v){
	mpc_add_si(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const unsigned int v){
	mpc_add_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const unsigned long v){
	mpc_add_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline const mpcomplex operator+(const mpcomplex& a, const mpcomplex& b){
	mpcomplex c;
	mpc_add(c.value, a.value, b.value, mpcomplex::default_rnd_mode);
	return c;
}

inline mpcomplex& mpcomplex::operator-=(const mpcomplex& v){
	mpc_sub(value, value, v.value, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const mpc_t v){
	mpc_sub(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const int v){
	mpc_add_si(value, value, -v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const long v){
	mpc_add_si(value, value, -v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const unsigned int v){
	mpc_sub_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const unsigned long v){
	mpc_sub_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline const mpcomplex mpcomplex::operator-()const{
	mpcomplex u(*this);
	mpc_neg(u.value, u.value, mpcomplex::default_rnd_mode);
	return u;
}

inline const mpcomplex operator-(const mpcomplex& a, const mpcomplex& b){
	mpcomplex c;
	mpc_sub(c.value, a.value, b.value, mpcomplex::default_rnd_mode);
	return c;
}

inline mpcomplex& mpcomplex::operator*=(const mpcomplex& v){
	mpc_mul(value, value, v.value, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const mpc_t v){
	mpc_mul(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const int v){
	mpc_mul_si(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const long v){
	mpc_mul_si(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const unsigned int v){
	mpc_mul_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const unsigned long v){
	mpc_mul_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline const mpcomplex operator*(const mpcomplex& a, const mpcomplex& b){
	mpcomplex c;
	mpc_mul(c.value, a.value, b.value, mpcomplex::default_rnd_mode);
	return c;
}

inline mpcomplex& mpcomplex::operator/=(const mpcomplex& v){
	mpc_div(value, value, v.value, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const mpc_t v){
	mpc_div(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const int v){
	mpc_div_ui(value, value, std::abs(v), mpcomplex::default_rnd_mode);
	if(v < 0) mpc_neg(value, value, rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const long v){
	mpc_div_ui(value, value, std::abs(v), mpcomplex::default_rnd_mode);
	if(v < 0) mpc_neg(value, value, rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const unsigned int v){
	mpc_div_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const unsigned long v){
	mpc_div_ui(value, value, v, mpcomplex::default_rnd_mode);
	return *this;
}

inline const mpcomplex operator/(const mpcomplex& a, const mpcomplex& b){
	mpcomplex c;
	mpc_div(c.value, a.value, b.value, mpcomplex::default_rnd_mode);
	return c;
}

inline mpcomplex& mpcomplex::operator<<=(const int u){
	mpc_mul_2si(value, value, u, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const long u){
	mpc_mul_2si(value, value, u, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const unsigned int u){
	mpc_mul_2ui(value, value, u, mpcomplex::default_rnd_mode);
	return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const unsigned long u){
	mpc_mul_2ui(value, value, u, mpcomplex::default_rnd_mode);
	return *this;
}

inline bool operator == (const mpcomplex& a, const mpcomplex& b){return (mpc_cmp(a.value,b.value)==0);}
inline bool operator != (const mpcomplex& a, const mpcomplex& b){return !(mpc_cmp(a.value,b.value)==0);}

inline std::string mpcomplex::to_string(int digits, int base) const{
	char* cstr = mpc_get_str(base, std::max(digits,0), value, MPC_RNDNN);
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
			if(mpfHalf < tolerance) halves[i-1].clear();
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
						halves[i-1].insert(exp+2, ".");
					} else {
						halves[i-1].erase(1, 1);
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

inline mpfr::mpreal mpcomplex::realPart(){
	mpfr::mpreal realPart;
	mpc_real(realPart.mpfr_ptr(), value, (mpfr_rnd_t)this->rnd_mode);
	return realPart;
}
inline mpfr::mpreal mpcomplex::imPart(){
	mpfr::mpreal imPart;
	mpc_imag(imPart.mpfr_ptr(), value, (mpfr_rnd_t)this->rnd_mode);
	return imPart;
}

inline bool mpcomplex::isReal(){
	if(this->imPart() < tolerance) return true;
	return false;
}

inline mpcomplex mpcomplex::sqr(){
	mpcomplex o(*this);
	mpc_sqr(o.value, o.value, rnd_mode);
	return o;
}

inline mpcomplex mpcomplex::sqrt(){
	mpcomplex o(*this);
	mpc_sqrt(o.value, o.value, rnd_mode);
	return o;
}

/*inline mpfr_t mpcomplex::abs(){
	mpfr_t o;
	mpc_abs(o, value, rnd_mode);
	return o;
}*/

inline const mpcomplex sqrt(const mpcomplex& x){
	mpcomplex o;
	mpc_sqrt(o.value, x.value, x.rnd_mode);
	return o;
}

inline const mpfr::mpreal abs(const mpcomplex& v){
	mpfr::mpreal abs;
	mpc_abs(abs.mpfr_ptr(), v.value, (mpfr_rnd_t)v.rnd_mode);
	return abs;
}

inline const std::string to_string(const mpcomplex& v, int digits = 0){
	return v.to_string(digits);
}

inline void swap(mpcomplex& a, mpcomplex& b){	mpc_swap(a.value, b.value);	}

}
#endif
