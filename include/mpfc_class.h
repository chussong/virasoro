#ifndef MPFC_CLASS_H_
#define MPFC_CLASS_H_

#include <string>
#include <mpc.h>

class mpfc_class{

	public:
		static mpc_rnd_t default_rnd_mode;
		static mpfr_prec_t default_prec;

		static void set_default_prec(const mpfr_prec_t prec);
		static void set_default_rnd_mode(const mpc_rnd_t rnd_mode);

		mpc_rnd_t rnd_mode;
		mpc_t value;

		mpfc_class();
		mpfc_class(const mpfc_class& u);
		mpfc_class(mpfc_class&& u);
		mpfc_class(long rlPart, long imPart = 0, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		mpfc_class(mpf_class rlPart, mpf_class imPart = 0, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		mpfc_class(std::string& str, int base = 10, mpfr_prec_t prec = default_prec, mpc_rnd_t rnd_mode = default_rnd_mode);
		~mpfc_class();
		
		mpfc_class& operator=(const mpfc_class& v);
		mpfc_class& operator=(mpfc_class&& v);
		mpfc_class& operator=(const mpc_t v);
		mpfc_class& operator=(const int v);
		mpfc_class& operator=(const long v);
		mpfc_class& operator=(const unsigned int v);
		mpfc_class& operator=(const unsigned long v);
		mpfc_class& operator=(const char* s);
		mpfc_class& operator=(const std::string& s);

		mpfc_class& operator+=(const mpfc_class& v);
		mpfc_class& operator+=(const mpc_t v);
		mpfc_class& operator+=(const int v);
		mpfc_class& operator+=(const long v);
		mpfc_class& operator+=(const unsigned int v);
		mpfc_class& operator+=(const unsigned long v);

		mpfc_class& operator-=(const mpfc_class& v);
		mpfc_class& operator-=(const mpc_t v);
		mpfc_class& operator-=(const int v);
		mpfc_class& operator-=(const long v);
		mpfc_class& operator-=(const unsigned int v);
		mpfc_class& operator-=(const unsigned long v);
		const mpfc_class operator-() const;

		mpfc_class& operator*=(const mpfc_class& v);
		mpfc_class& operator*=(const mpc_t v);
		mpfc_class& operator*=(const int v);
		mpfc_class& operator*=(const long v);
		mpfc_class& operator*=(const unsigned int v);
		mpfc_class& operator*=(const unsigned long v);

		mpfc_class& operator/=(const mpfc_class& v);
		mpfc_class& operator/=(const mpc_t v);
		mpfc_class& operator/=(const int v);
		mpfc_class& operator/=(const long v);
		mpfc_class& operator/=(const unsigned int v);
		mpfc_class& operator/=(const unsigned long v);

		// <<= is fast multiplication by 2^u
		mpfc_class& operator<<=(const int u);
		mpfc_class& operator<<=(const long u);
		mpfc_class& operator<<=(const unsigned int u);
		mpfc_class& operator<<=(const unsigned long u);

		std::string to_string(int digits = 0, int base = 10) const;

		mpfc_class sqr();
		mpfc_class sqrt();

		const mpfc_class sqrt(const mpfc_class& v);
		const mpfr_t* abs(const mpfc_class& v);

		friend void swap(mpfc_class& x, mpfc_class& y);
};

mpc_rnd_t mpfc_class::default_rnd_mode = MPC_RNDZZ;
mpfr_prec_t mpfc_class::default_prec = 64;

inline void mpfc_class::set_default_prec(const mpfr_prec_t prec){
	default_prec = prec;
}

inline void mpfc_class::set_default_rnd_mode(const mpc_rnd_t rnd_mode){
	default_rnd_mode = rnd_mode;
}

inline mpfc_class::mpfc_class(){
	rnd_mode = default_rnd_mode;
	mpc_init2(value, mpfc_class::default_prec);
	mpc_set_ui(value, 0, rnd_mode);
}

inline mpfc_class::mpfc_class(const mpfc_class& u){
	rnd_mode = u.rnd_mode;
	mpc_init2(value, mpfc_class::default_prec);
	mpc_set(value, u.value, rnd_mode);
}

inline mpfc_class::mpfc_class(mpfc_class&& u){
	rnd_mode = u.rnd_mode;
//	mpc_init2(value, mpfc_class::default_prec);		not sure if we need this or not
	mpc_swap(value, u.value);
}

inline mpfc_class::mpfc_class(long rlPart, long imPart, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_si_si(value, rlPart, imPart, rnd_mode);
}

inline mpfc_class::mpfc_class(mpf_class rlPart, mpf_class imPart, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_f_f(value, rlPart.get_mpf_t(), imPart.get_mpf_t(), rnd_mode);
}

inline mpfc_class::mpfc_class(std::string& str, int base, mpfr_prec_t prec, mpc_rnd_t rnd_mode):rnd_mode(rnd_mode){
	mpc_init2(value, prec);
	mpc_set_str(value, str.c_str(), base, rnd_mode);
}

inline mpfc_class::~mpfc_class(){
	mpc_clear(value);
}

// this is needed for the templates
namespace internal{

	template <typename ArgumentType> struct result_type {};

	template <> struct result_type<mpfc_class>			{typedef mpfc_class type;};
	template <> struct result_type<int>					{typedef mpfc_class type;};
	template <> struct result_type<long>				{typedef mpfc_class type;};
	template <> struct result_type<unsigned int>		{typedef mpfc_class type;};
	template <> struct result_type<unsigned long>		{typedef mpfc_class type;};
}

template <typename Rhs>
inline const typename internal::result_type<Rhs>::type
	operator+(const mpfc_class& lhs, const Rhs& rhs){ return mpfc_class(lhs) += rhs; }

template <typename Lhs>
inline const typename internal::result_type<Lhs>::type
	operator+(const Lhs& lhs, const mpfc_class& rhs){ return mpfc_class(rhs) += lhs; }

template <typename Rhs>
inline const typename internal::result_type<Rhs>::type
	operator-(const mpfc_class& lhs, const Rhs& rhs){ return mpfc_class(lhs) -= rhs; }

template <typename Lhs>
inline const typename internal::result_type<Lhs>::type
	operator-(const Lhs& lhs, const mpfc_class& rhs){ return mpfc_class(lhs) -= rhs; }

template <typename Rhs>
inline const typename internal::result_type<Rhs>::type
	operator*(const mpfc_class& lhs, const Rhs& rhs){ return mpfc_class(lhs) *= rhs; }

template <typename Lhs>
inline const typename internal::result_type<Lhs>::type
	operator*(const Lhs& lhs, const mpfc_class& rhs){ return mpfc_class(rhs) *= lhs; }

template <typename Rhs>
inline const typename internal::result_type<Rhs>::type
	operator/(const mpfc_class& lhs, const Rhs& rhs){ return mpfc_class(lhs) /= rhs; }

template <typename Lhs>
inline const typename internal::result_type<Lhs>::type
	operator/(const Lhs& lhs, const mpfc_class& rhs){ return mpfc_class(lhs) /= rhs; }

inline mpfc_class& mpfc_class::operator=(const mpfc_class& v){
	mpc_set(value, v.value, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(mpfc_class&& v){
	mpc_swap(value, v.value);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const mpc_t v){
	mpc_set(value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const int v){
	mpc_set_si(value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const long v){
	mpc_set_si(value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const unsigned int v){
	mpc_set_ui(value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const unsigned long v){
	mpc_set_ui(value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const char* s){
	mpc_t t;
	mpc_init2(t, mpfc_class::default_rnd_mode);

	if(0 == mpc_set_str(t, s, 10, mpfc_class::default_rnd_mode)){
		mpc_set(value, t, mpfc_class::default_rnd_mode);
	}
	mpc_clear(t);
	return *this;
}

inline mpfc_class& mpfc_class::operator=(const std::string& s){
	mpc_t t;
	mpc_init2(t, mpfc_class::default_rnd_mode);

	if(0 == mpc_set_str(t, s.c_str(), 10, mpfc_class::default_rnd_mode)){
		mpc_set(value, t, mpfc_class::default_rnd_mode);
	}
	mpc_clear(t);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const mpfc_class& v){
	mpc_add(value, value, v.value, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const mpc_t v){
	mpc_add(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const int v){
	mpc_add_si(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const long v){
	mpc_add_si(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const unsigned int v){
	mpc_add_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator+=(const unsigned long v){
	mpc_add_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline const mpfc_class operator+(const mpfc_class& a, const mpfc_class& b){
	mpfc_class c;
	mpc_add(c.value, a.value, b.value, mpfc_class::default_rnd_mode);
	return c;
}

inline mpfc_class& mpfc_class::operator-=(const mpfc_class& v){
	mpc_sub(value, value, v.value, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator-=(const mpc_t v){
	mpc_sub(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator-=(const int v){
	mpc_add_si(value, value, -v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator-=(const long v){
	mpc_add_si(value, value, -v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator-=(const unsigned int v){
	mpc_sub_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator-=(const unsigned long v){
	mpc_sub_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline const mpfc_class mpfc_class::operator-()const{
	mpfc_class u(*this);
	mpc_neg(u.value, u.value, mpfc_class::default_rnd_mode);
	return u;
}

inline const mpfc_class operator-(const mpfc_class& a, const mpfc_class& b){
	mpfc_class c;
	mpc_sub(c.value, a.value, b.value, mpfc_class::default_rnd_mode);
	return c;
}

inline mpfc_class& mpfc_class::operator*=(const mpfc_class& v){
	mpc_mul(value, value, v.value, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator*=(const mpc_t v){
	mpc_mul(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator*=(const int v){
	mpc_mul_si(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator*=(const long v){
	mpc_mul_si(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator*=(const unsigned int v){
	mpc_mul_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator*=(const unsigned long v){
	mpc_mul_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline const mpfc_class operator*(const mpfc_class& a, const mpfc_class& b){
	mpfc_class c;
	mpc_mul(c.value, a.value, b.value, mpfc_class::default_rnd_mode);
	return c;
}

inline mpfc_class& mpfc_class::operator/=(const mpfc_class& v){
	mpc_div(value, value, v.value, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator/=(const mpc_t v){
	mpc_div(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator/=(const int v){
	mpc_div_ui(value, value, std::abs(v), mpfc_class::default_rnd_mode);
	if(v < 0) mpc_neg(value, value, rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator/=(const long v){
	mpc_div_ui(value, value, std::abs(v), mpfc_class::default_rnd_mode);
	if(v < 0) mpc_neg(value, value, rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator/=(const unsigned int v){
	mpc_div_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator/=(const unsigned long v){
	mpc_div_ui(value, value, v, mpfc_class::default_rnd_mode);
	return *this;
}

inline const mpfc_class operator/(const mpfc_class& a, const mpfc_class& b){
	mpfc_class c;
	mpc_div(c.value, a.value, b.value, mpfc_class::default_rnd_mode);
	return c;
}

inline mpfc_class& mpfc_class::operator<<=(const int u){
	mpc_mul_2si(value, value, u, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator<<=(const long u){
	mpc_mul_2si(value, value, u, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator<<=(const unsigned int u){
	mpc_mul_2ui(value, value, u, mpfc_class::default_rnd_mode);
	return *this;
}

inline mpfc_class& mpfc_class::operator<<=(const unsigned long u){
	mpc_mul_2ui(value, value, u, mpfc_class::default_rnd_mode);
	return *this;
}

inline bool operator == (const mpfc_class& a, const mpfc_class& b){return (mpc_cmp(a.value,b.value)==0);}
inline bool operator != (const mpfc_class& a, const mpfc_class& b){return !(mpc_cmp(a.value,b.value)==0);}

inline std::string mpfc_class::to_string(int digits, int base) const{
	char* cstr = mpc_get_str(base, digits, value, rnd_mode);
	std::string output = cstr;
	mpc_free_str(cstr);
	return output;
}

inline mpfc_class mpfc_class::sqr(){
	mpfc_class o(*this);
	mpc_sqr(o.value, o.value, rnd_mode);
	return o;
}

inline mpfc_class mpfc_class::sqrt(){
	mpfc_class o(*this);
	mpc_sqrt(o.value, o.value, rnd_mode);
	return o;
}

/*inline mpfr_t mpfc_class::abs(){
	mpfr_t o;
	mpc_abs(o, value, rnd_mode);
	return o;
}*/

inline const mpfc_class sqrt(const mpfc_class& x){
	mpfc_class o;
	mpc_sqrt(o.value, x.value, x.rnd_mode);
	return o;
}

inline const mpf_class abs(const mpfc_class& v){
	mpfr_t abs;
	mpfr_init(abs);
	mpc_abs(abs, v.value, (mpfr_rnd_t)v.rnd_mode);
	mpf_class ret;
	mpfr_get_f(ret.get_mpf_t(), abs, (mpfr_rnd_t)v.rnd_mode);
	return ret;
}

inline void swap(mpfc_class& a, mpfc_class& b){	mpc_swap(a.value, b.value);	}

#endif
