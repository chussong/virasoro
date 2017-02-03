/*
    MPFR C++: Multi-precision floating point number class for C++. 
    Based on MPFR library:    http://mpfr.org

    Project homepage:    http://www.holoborodko.com/pavel/mpfr
    Contact e-mail:      pavel@holoborodko.com

    Copyright (c) 2008-2015 Pavel Holoborodko

    Contributors:
    Dmitriy Gubanov, Konstantin Holoborodko, Brian Gladman,
    Charles Hussong, Helmut Jarausch, Fokko Beekhof, Ulrich Mutze,
    Heinz van Saanen, Pere Constans, Peter van Hoof, Gael Guennebaud,
    Tsai Chia Cheng, Alexei Zubanov, Jauhien Piatlicki, Victor Berger,
    John Westwood, Petr Aleksandrov, Orion Poplawski, Charles Karney,
    Arash Partow, Rodney James, Jorge Leitao.

    Licensing:
    (A) MPFR C++ is under GNU General Public License ("GPL").
    
    (B) Non-free licenses may also be purchased from the author, for users who 
        do not want their programs protected by the GPL.

        The non-free licenses are for users that wish to use MPFR C++ in 
        their products but are unwilling to release their software 
        under the GPL (which would require them to release source code 
        and allow free redistribution).

        Such users can purchase an unlimited-use license from the author.
        Contact us for more details.
    
    GNU General Public License ("GPL") copyright permissions statement:
    **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __MPCOMPLEX_H__
#define __MPCOMPLEX_H__

#include "mpreal.h"
#include <mpc.h>

#ifdef MPREAL_HAVE_MOVE_SUPPORT
    #define mpc_is_initialized(x)      (0 != (x)->re->_mpfr_d || 0 != (x)->im->_mpfr_d)
#endif

#define defround (mpc_rnd_t)MPC_RND(mpreal::get_default_rnd(), mpreal::get_default_rnd())

// forward declarations in order to put math functions in mpfr namespace where mpreal is
namespace mpfr{
	using std::complex;
	const mpreal          abs  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
	const mpreal          real (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
	const mpreal          imag (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);

    const complex<mpreal> exp  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround); 
    const complex<mpreal> log  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> log10(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);

    const complex<mpreal> pow (const complex<mpreal>& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const mpreal& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const complex<mpreal>& a, const mpreal& b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const complex<mpreal>& a, const unsigned long int b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const complex<mpreal>& a, const unsigned int b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const complex<mpreal>& a, const long int b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const complex<mpreal>& a, const int b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const unsigned long int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> pow (const unsigned long int a, const unsigned long int b, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqr (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrtcomp(const mpreal& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const unsigned long int v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const unsigned int v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const long int v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const int v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const long double v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sqrt(const double v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> root(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode = defround);

    inline const complex<mpreal> mul_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode = defround);
    inline const complex<mpreal> mul_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode = defround);
    inline const complex<mpreal> div_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode = defround);
    inline const complex<mpreal> div_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode = defround);
    
    const complex<mpreal> sin(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> cos(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> tan(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sec(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> csc(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> cot(const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    int sin_cos(complex<mpreal>& s, complex<mpreal>& c, const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);

    const complex<mpreal> asin  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acos  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> atan  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acot  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> asec  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acsc  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);

    const complex<mpreal> sinh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> cosh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> tanh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> sech  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> csch  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> coth  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> asinh (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acosh (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> atanh (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acoth (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> asech (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);
    const complex<mpreal> acsch (const complex<mpreal>& v, mpc_rnd_t rnd_mode = defround);

    const complex<mpreal> mpc_urandom (gmp_randstate_t& state);     // use gmp_randinit_default() to init state, gmp_randclear() to clear

    // constant infinity value in case someone cares
    const complex<mpreal> mpc_const_infinity(int ReSign, int ImSign, mp_prec_t prec);
}
namespace std{
using mpfr::mpreal;

template<>
class complex<mpreal> {
private:
    mpc_t mp;
    
public:
    
    // MPC doesn't have default rnd/prec, so we use mpreal's defaults for both components
    inline static mpc_rnd_t   get_default_rnd()   {    return (mpc_rnd_t)MPC_RND(mpreal::get_default_rnd(), mpreal::get_default_rnd());       }
    inline static mpfr_prec_t  get_default_prec() {    return mpreal::get_default_prec(); }

    // Constructors && type conversions
    complex<mpreal>();
    complex<mpreal>(const complex<mpreal>& u);
    complex<mpreal>(const complex<mpreal>& u, mpc_rnd_t mode);
    // Construct complex<mpreal> from mpc_t structure.
    // shared = true allows to avoid deep copy, so that complex<mpreal> and 'u' share the same data & pointers.    
    complex<mpreal>(const mpc_t& u, bool shared = false);   

    // these are templates so that we can try constructing through mpreal if no mpc_set_ exists
    template<class X> explicit complex<mpreal>(const X re, mp_prec_t prec = get_default_prec(), mpc_rnd_t rnd_mode = get_default_rnd());
    template<class X> complex<mpreal>(const X re, const X im, mp_prec_t prec = get_default_prec(), mpc_rnd_t rnd_mode = get_default_rnd());
    template<class X> explicit complex<mpreal>(const std::complex<X> u, mp_prec_t prec = get_default_prec(), mpc_rnd_t mode = get_default_rnd());

#ifdef _COMPLEX_H
    explicit complex<mpreal>(const double _Complex u, mp_prec_t prec = get_default_prec(), mpc_rnd_t mode = get_default_rnd());
    explicit complex<mpreal>(const long double _Complex u, mp_prec_t prec = get_default_prec(), mpc_rnd_t mode = get_default_rnd());
    complex<mpreal>& operator=(const double _Complex v);
    complex<mpreal>& operator=(const long double _Complex v);

    _Complex float              toCFloat    (mpc_rnd_t mode = MPC_RNDNN)    const;
    _Complex double             toCDouble   (mpc_rnd_t mode = MPC_RNDNN)    const;
    _Complex long double        toCLDouble  (mpc_rnd_t mode = MPC_RNDNN)    const;
    explicit operator _Complex float             () const { return toCFloat();               }
    explicit operator _Complex double            () const { return toCDouble();              }
    explicit operator _Complex long double       () const { return toCLDouble();             }
#endif
    
    explicit complex<mpreal>(const char* s,             mp_prec_t prec = complex<mpreal>::get_default_prec(), int base = 10, mpc_rnd_t mode = complex<mpreal>::get_default_rnd());
    explicit complex<mpreal>(const std::string& s,      mp_prec_t prec = complex<mpreal>::get_default_prec(), int base = 10, mpc_rnd_t mode = complex<mpreal>::get_default_rnd());

    ~complex<mpreal>();                           

#ifdef MPREAL_HAVE_MOVE_SUPPORT
    complex<mpreal>& operator=(complex<mpreal>&& v);
    complex<mpreal>(complex<mpreal>&& u);
#endif

    // Operations
    // =
    // +, -, *, /, ++, --, <<, >> 
    // *=, +=, -=, /=,
    // <, >, ==, <=, >=

    // =
    complex<mpreal>& operator=(const complex<mpreal>& v);
    template<typename real_t> complex<mpreal>& operator= (const std::complex<real_t>& z);
    template<class X> complex<mpreal>& operator=(const X v);
    template<class X> complex<mpreal>& real(const X v);
    template<class X> complex<mpreal>& imag(const X v);
    complex<mpreal>& operator=(const char* s);
    complex<mpreal>& real(const char* s);
    complex<mpreal>& imag(const char* s);
    complex<mpreal>& operator=(const std::string& s);
    complex<mpreal>& real(const std::string& s);
    complex<mpreal>& imag(const std::string& s);

    complex<mpreal>& operator+=(const complex<mpreal>& v);
    complex<mpreal>& operator+=(const mpreal &v);
    complex<mpreal>& operator+=(const mpf_t v);
    complex<mpreal>& operator+=(const long double u);
    complex<mpreal>& operator+=(const double u);
    complex<mpreal>& operator+=(const unsigned long int u);
    complex<mpreal>& operator+=(const unsigned int u);
    complex<mpreal>& operator+=(const long int u);
    complex<mpreal>& operator+=(const int u);

    complex<mpreal>& operator+=(const long long int  u);
    complex<mpreal>& operator+=(const unsigned long long int u);
    complex<mpreal>& operator-=(const long long int  u);
    complex<mpreal>& operator-=(const unsigned long long int u);
    complex<mpreal>& operator*=(const long long int  u);
    complex<mpreal>& operator*=(const unsigned long long int u);
    complex<mpreal>& operator/=(const long long int  u);
    complex<mpreal>& operator/=(const unsigned long long int u);

    const complex<mpreal> operator+() const;
    complex<mpreal>& operator++ ();
    const complex<mpreal>  operator++ (int); 

    complex<mpreal>& operator-=(const complex<mpreal>& v);
    complex<mpreal>& operator-=(const mpreal& v);
    complex<mpreal>& operator-=(const long double u);
    complex<mpreal>& operator-=(const double u);
    complex<mpreal>& operator-=(const unsigned long int u);
    complex<mpreal>& operator-=(const unsigned int u);
    complex<mpreal>& operator-=(const long int u);
    complex<mpreal>& operator-=(const int u);
    const complex<mpreal> operator-() const;
    friend complex<mpreal> operator-(const mpreal b           ,  const complex<mpreal>& a);
    friend complex<mpreal> operator-(const unsigned long int b,  const complex<mpreal>& a);
    friend complex<mpreal> operator-(const unsigned int b,       const complex<mpreal>& a);
    friend complex<mpreal> operator-(const unsigned short int b, const complex<mpreal>& a);
    complex<mpreal>& operator-- ();    
    const complex<mpreal>  operator-- (int);

    complex<mpreal>& operator*=(const complex<mpreal>& v);
    complex<mpreal>& operator*=(const mpreal& v);
    complex<mpreal>& operator*=(const long double v);
    complex<mpreal>& operator*=(const double v);
    complex<mpreal>& operator*=(const unsigned long int v);
    complex<mpreal>& operator*=(const unsigned int v);
    complex<mpreal>& operator*=(const long int v);
    complex<mpreal>& operator*=(const int v);
    
    complex<mpreal>& operator/=(const complex<mpreal>& v);
    complex<mpreal>& operator/=(const mpreal& v);
    complex<mpreal>& operator/=(const long double v);
    complex<mpreal>& operator/=(const double v);
    complex<mpreal>& operator/=(const unsigned long int v);
    complex<mpreal>& operator/=(const unsigned int v);
    complex<mpreal>& operator/=(const long int v);
    complex<mpreal>& operator/=(const int v);
    friend complex<mpreal> operator/(const mpreal b           ,  const complex<mpreal>& a);
    friend complex<mpreal> operator/(const unsigned long int b,  const complex<mpreal>& a);
    friend complex<mpreal> operator/(const unsigned int b,       const complex<mpreal>& a);
    friend complex<mpreal> operator/(const unsigned short int b, const complex<mpreal>& a);

    //<<= Fast Multiplication by 2^u
    complex<mpreal>& operator<<=(const unsigned long int u);
    complex<mpreal>& operator<<=(const unsigned int u);
    complex<mpreal>& operator<<=(const long int u);
    complex<mpreal>& operator<<=(const int u);

    //>>= Fast Division by 2^u
    complex<mpreal>& operator>>=(const unsigned long int u);
    complex<mpreal>& operator>>=(const unsigned int u);
    complex<mpreal>& operator>>=(const long int u);
    complex<mpreal>& operator>>=(const int u);

    // Type Conversion operators
    std::complex<float>         toFloat     (mpc_rnd_t mode = MPC_RNDNN)    const;
    std::complex<double>        toDouble    (mpc_rnd_t mode = MPC_RNDNN)    const;
    std::complex<long double>   toLDouble   (mpc_rnd_t mode = MPC_RNDNN)    const;

#if defined (MPREAL_HAVE_EXPLICIT_CONVERTERS)
    explicit operator std::complex<float>        () const { return toFloat();                }
    explicit operator std::complex<double>       () const { return toDouble();               }
    explicit operator std::complex<long double>  () const { return toLDouble();              }
#endif

    // Get raw pointers so that complex<mpreal> can be directly used in raw mpc_* functions
    ::mpc_ptr     mpc_ptr();
    ::mpc_srcptr  mpc_ptr()       const;
    ::mpc_srcptr  mpc_srcptr()    const;
    ::mpfr_ptr    mpc_re_ptr()          { return mpc_realref(mpc_ptr());    }
    ::mpfr_srcptr mpc_re_ptr()    const { return mpc_realref(mpc_srcptr()); }
    ::mpfr_srcptr mpc_re_srcptr() const { return mpc_realref(mpc_srcptr()); }
    ::mpfr_ptr    mpc_im_ptr()          { return mpc_imagref(mpc_ptr());    }
    ::mpfr_srcptr mpc_im_ptr()    const { return mpc_imagref(mpc_srcptr()); }
    ::mpfr_srcptr mpc_im_srcptr() const { return mpc_imagref(mpc_srcptr()); }

    // Convert mpcomplex to string with n significant digits in base b
    // n = -1 -> convert with the maximum available digits
    std::string toString(int n = -1, int b = 10, mpc_rnd_t mode = complex<mpreal>::get_default_rnd()) const;

#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
    std::string mpfrString(const ::mpfr_srcptr n, const std::string& format) const;
#endif

    std::ostream& output(std::ostream& os) const;

    mpreal real() const;
    mpreal imag() const;

	// std::complex function specializations
/*    friend const mpreal abs  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal norm (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal real (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal imag (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal arg  (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());

    friend const complex<mpreal> conj (const complex<mpreal>& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const complex<mpreal> proj (const complex<mpreal>& z, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const complex<mpreal> polar(const mpreal& r, const mpreal& theta = 0, mpc_rnd_t rnd_mode = get_default_rnd());

    friend const complex<mpreal> exp  (const complex<mpreal>& v, mpc_rnd_t rnd_mode); 
    friend const complex<mpreal> log  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> log10(const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> pow (const complex<mpreal>& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> pow (const mpreal& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> pow (const complex<mpreal>& a, const mpreal& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> sqrt(const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> sin(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> cos(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> tan(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> asin  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> acos  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> atan  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> sinh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> cosh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> tanh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> asinh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> acosh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> atanh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
*/
    // Math Functions
	friend const mpreal          mpfr::abs  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
	friend const mpreal          mpfr::real (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
	friend const mpreal          mpfr::imag (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::exp  (const complex<mpreal>& v, mpc_rnd_t rnd_mode); 
    friend const complex<mpreal> mpfr::log  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::log10(const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const mpreal& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const mpreal& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const unsigned long int b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const unsigned int b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const long int b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const complex<mpreal>& a, const int b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const unsigned long int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::pow (const unsigned long int a, const unsigned long int b, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqr (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrtcomp(const mpreal& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const unsigned long int v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const unsigned int v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const long int v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const int v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const long double v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sqrt(const double v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> mpfr::mul_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::mul_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::div_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::div_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode);
    
    friend const complex<mpreal> mpfr::sin(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::cos(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::tan(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sec(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::csc(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::cot(const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend int mpfr::sin_cos(complex<mpreal>& s, complex<mpreal>& c, const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> mpfr::asin  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acos  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::atan  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acot  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::asec  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acsc  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> mpfr::sinh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::cosh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::tanh  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::sech  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::csch  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::coth  (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::asinh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acosh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::atanh (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acoth (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::asech (const complex<mpreal>& v, mpc_rnd_t rnd_mode);
    friend const complex<mpreal> mpfr::acsch (const complex<mpreal>& v, mpc_rnd_t rnd_mode);

    friend const complex<mpreal> mpfr::mpc_urandom (gmp_randstate_t& state);     // use gmp_randinit_default() to init state, gmp_randclear() to clear

    // constant infinity value in case someone cares
    friend const complex<mpreal> mpfr::mpc_const_infinity(int ReSign, int ImSign, mp_prec_t prec);

    // Output/ Input
    friend std::ostream& operator<<(std::ostream& os, const complex<mpreal>& v);
    friend std::istream& operator>>(std::istream& is, complex<mpreal>& v);

    // Set/Get instance properties
    inline mp_prec_t    get_prec() const;
    inline void         set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode = get_default_rnd());    // Change precision with rounding mode
    inline int          getPrecision() const;
    inline complex<mpreal>&   setPrecision(int Precision, mpc_rnd_t RoundingMode = get_default_rnd());

    // Set mpcomplex to +/- inf, NaN, +/-0
    complex<mpreal>&        setInf  (int ReSign = +1, int ImSign = 0);    
    complex<mpreal>&        setNan  ();
    complex<mpreal>&        setZero (int ReSign = +1, int ImSign = +1);
    complex<mpreal>&        setSigns(int ReSign, int ImSign, mpc_rnd_t rnd_mode);

    // Set/Get global properties
    static void            set_default_rnd(mpc_rnd_t rnd_mode);

    static mp_exp_t  get_emin (void);
    static mp_exp_t  get_emax (void);
    static mp_exp_t  get_emin_min (void);
    static mp_exp_t  get_emin_max (void);
    static mp_exp_t  get_emax_min (void);
    static mp_exp_t  get_emax_max (void);
    static int       set_emin (mp_exp_t exp);
    static int       set_emax (mp_exp_t exp);

    // Efficient swapping of two complex<mpreal> values - needed for std algorithms
    friend void swap(complex<mpreal>& x, complex<mpreal>& y);
    
private:
    // Human friendly Debug Preview in Visual Studio.
    // Put one of these lines:
    //
    // std::complex<mpfr::mpreal>=<DebugView>                              ; Show value only
    // std::complex<mpfr::mpreal>=<DebugView>, <mp[0]._mpfr_prec,u>bits    ; Show value & precision
    // 
    // at the beginning of
    // [Visual Studio Installation Folder]\Common7\Packages\Debugger\autoexp.dat
    MPREAL_MSVC_DEBUGVIEW_DATA

    // "Smart" resources deallocation. Checks if instance initialized before deletion.
    void clear(::mpc_ptr);
};

//////////////////////////////////////////////////////////////////////////
// Constructors & converters
// Default constructor: creates mp number and initializes it to 0.

inline complex<mpreal>::complex() 
{ 
    mpc_init2(mpc_ptr(), complex<mpreal>::get_default_prec()); 
    mpc_set_ui(mpc_ptr(), 0, complex<mpreal>::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline complex<mpreal>::complex(const complex<mpreal>& u) 
{
    mpc_init2(mpc_ptr(), mpc_get_prec(u.mpc_srcptr()));
    mpc_set  (mpc_ptr(), u.mpc_srcptr(), complex<mpreal>::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline complex<mpreal>::complex(const complex<mpreal>& u, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), mpc_get_prec(u.mpc_srcptr()));
    mpc_set  (mpc_ptr(), u.mpc_srcptr(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

#ifdef MPREAL_HAVE_MOVE_SUPPORT
inline complex<mpreal>::complex(complex<mpreal>&& other)
{
//    mpc_set_uninitialized(mpc_ptr());     
	mpc_re_ptr()->_mpfr_d = 0;      // make sure "other" holds no pointer to actual data 
	mpc_im_ptr()->_mpfr_d = 0;
    mpc_swap(mpc_ptr(), other.mpc_ptr());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline complex<mpreal>& complex<mpreal>::operator=(complex<mpreal>&& other)
{
    mpc_swap(mpc_ptr(), other.mpc_ptr());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}
#endif

inline complex<mpreal>::complex(const mpc_t& u, bool shared)
{
    if(shared)
    {
        std::memcpy(mpc_ptr(), u, sizeof(mpc_t));
    }
    else
    {
        mpc_init2(mpc_ptr(), mpc_get_prec(u));
        mpc_set  (mpc_ptr(), u, get_default_rnd());
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

// if none of the subsequent specializations work, try to make an mpreal out of re
template<class X>
inline complex<mpreal>::complex(const X re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpreal re2(re);
    mpc_init2(mpc_ptr(), prec);
    mpc_set_fr(mpc_ptr(), re2.mpfr_srcptr(), rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpreal& re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), std::max(mpfr_get_prec(re.mpfr_srcptr()), static_cast<long int>(prec)));
    mpc_set_fr(mpc_ptr(), re.mpfr_srcptr(), rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned long int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned long long int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_uj(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long long int re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_sj(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const double re, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);

#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(re, MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpc_set_d(mpc_ptr(), re, mode);
    }else
        throw conversion_overflow();
#else
    mpc_set_d(mpc_ptr(), re, mode);
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long double re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ld(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpz_t& re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_z(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpq_t& re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_q(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpf_t& re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), std::max(static_cast<long int>(mpf_get_prec(re)), prec));
    mpc_set_f(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpfr_t& re, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_fr(mpc_ptr(), re, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

// if none of the subsequent specializations work, try to make mpreals out of re and im
template<class X>
inline complex<mpreal>::complex(const X re, const X im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpreal re2(re);
    mpreal im2(im);
    mpc_init2(mpc_ptr(), prec);
    mpc_set_fr_fr(mpc_ptr(), re2.mpfr_srcptr(), im2.mpfr_srcptr(), rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpreal& re, const mpreal& im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), std::max(std::max(mpfr_get_prec(re.mpfr_srcptr()), mpfr_get_prec(im.mpfr_srcptr())), static_cast<long int>(prec)));
    mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), im.mpfr_srcptr(), rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned int re, const unsigned int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui_ui(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const int re, const int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si_si(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned long int re, const unsigned long int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui_ui(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long int re, const long int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si_si(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const unsigned long long int re, const unsigned long long int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_uj_uj(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long long int re, const long long int im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_sj_sj(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const double re, const double im, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);

#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(re, MPREAL_DOUBLE_BITS_OVERFLOW) && mpfr::mpreal::fits_in_bits(im, MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpc_set_d_d(mpc_ptr(), re, im, mode);
    }else
        throw conversion_overflow();
#else
    mpc_set_d_d(mpc_ptr(), re, im, mode);
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long double re, const long double im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ld_ld(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpz_t& re, const mpz_t& im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_z_z(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpq_t& re, const mpq_t& im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_q_q(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpf_t& re, const mpf_t& im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), std::max(std::max(mpf_get_prec(re), mpf_get_prec(im)), static_cast<long unsigned int>(prec)));
    mpc_set_f_f(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const mpfr_t& re, const mpfr_t& im, mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_fr_fr(mpc_ptr(), re, im, rnd_mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

// if none of the subsequent specializations work, try to make mpreals out of u
template<class X>
inline complex<mpreal>::complex(const std::complex<X> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpreal re(u.real());
    mpreal im(u.imag());
    mpc_init2(mpc_ptr(), prec);
    mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), im.mpfr_srcptr(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<unsigned int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui_ui(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si_si(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<unsigned long int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ui_ui(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<long int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_si_si(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<unsigned long long int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_uj_uj(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<long long int> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_sj_sj(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<double> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);

#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(u.real(), MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpc_set_d_d(mpc_ptr(), u.real(), u.imag(), mode);
    }else
        throw conversion_overflow();
#else
    mpc_set_d_d(mpc_ptr(), u.real(), u.imag(), mode);
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const std::complex<long double> u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ld_ld(mpc_ptr(), u.real(), u.imag(), mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

// c-style complex numbers are a special case
#ifdef _COMPLEX_H
template<>
inline complex<mpreal>::complex(const double _Complex u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_dc(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

template<>
inline complex<mpreal>::complex(const long double _Complex u, mp_prec_t prec, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), prec);
    mpc_set_ldc(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}
#endif

inline complex<mpreal>::complex(const char* s, mp_prec_t prec, int base, mpc_rnd_t mode)
{
    mpc_init2  (mpc_ptr(), prec);
    mpc_set_str(mpc_ptr(), s, base, mode); 

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline complex<mpreal>::complex(const std::string& s, mp_prec_t prec, int base, mpc_rnd_t mode)
{
    mpc_init2  (mpc_ptr(), prec);
    mpc_set_str(mpc_ptr(), s.c_str(), base, mode); 

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline void complex<mpreal>::clear(::mpc_ptr x)
{
#ifdef MPREAL_HAVE_MOVE_SUPPORT
    if(mpc_is_initialized(x)) 
#endif
    mpc_clear(x);
}

inline complex<mpreal>::~complex<mpreal>() 
{ 
    clear(mpc_ptr());
}                           

// internal namespace needed for template magic
namespace complex_internal{

    // Use SFINAE to restrict arithmetic operations instantiation only for numeric types.
    // If operand is one of the below types, the template specializes to one containing
    // the type->mpcomplex map; if not, type does not exist and arithmetic template fails.
    // This is needed for smooth integration with libraries based on expression templates, like Eigen.
    // TODO: Do the same for boolean operators.
    template <typename ArgumentType> struct result_type {};    
    
    template <> struct result_type<double>              {typedef complex<mpreal> type;};    
    template <> struct result_type<unsigned long int>   {typedef complex<mpreal> type;};    
    template <> struct result_type<unsigned int>        {typedef complex<mpreal> type;};    
    template <> struct result_type<unsigned short int>  {typedef complex<mpreal> type;};    
    template <> struct result_type<long int>            {typedef complex<mpreal> type;};    
    template <> struct result_type<int>                 {typedef complex<mpreal> type;};    
    template <> struct result_type<short int>           {typedef complex<mpreal> type;};    
    template <> struct result_type<long long>           {typedef complex<mpreal> type;};    
    template <> struct result_type<unsigned long long>  {typedef complex<mpreal> type;};    
}

// + Addition
template <typename Rhs> 
inline typename complex_internal::result_type<Rhs>::type 
    operator+(const complex<mpreal>& lhs, const Rhs& rhs){ return complex<mpreal>(lhs) += rhs;    }

template <typename Lhs> 
inline typename complex_internal::result_type<Lhs>::type 
    operator+(const Lhs& lhs, const complex<mpreal>& rhs){ return complex<mpreal>(rhs) += lhs;    } 

// - Subtraction
template <typename Rhs> 
inline typename complex_internal::result_type<Rhs>::type 
    operator-(const complex<mpreal>& lhs, const Rhs& rhs){ return complex<mpreal>(lhs) -= rhs;    }

template <typename Lhs> 
inline typename complex_internal::result_type<Lhs>::type 
    operator-(const Lhs& lhs, const complex<mpreal>& rhs){ return complex<mpreal>(lhs) -= rhs;    }

// * Multiplication
template <typename Rhs> 
inline typename complex_internal::result_type<Rhs>::type 
    operator*(const complex<mpreal>& lhs, const Rhs& rhs){ return complex<mpreal>(lhs) *= rhs;    }

template <typename Lhs> 
inline typename complex_internal::result_type<Lhs>::type 
    operator*(const Lhs& lhs, const complex<mpreal>& rhs){ return complex<mpreal>(rhs) *= lhs;    } 

// / Division
template <typename Rhs> 
inline typename complex_internal::result_type<Rhs>::type 
    operator/(const complex<mpreal>& lhs, const Rhs& rhs){ return complex<mpreal>(lhs) /= rhs;    }

template <typename Lhs> 
inline typename complex_internal::result_type<Lhs>::type 
    operator/(const Lhs& lhs, const complex<mpreal>& rhs){ return complex<mpreal>(lhs) /= rhs;    }

/*} namespace mpfr {
// These functions are not part of the std::complex spec so they can't be in namespace std
//////////////////////////////////////////////////////////////////////////
// pow
const complex<mpreal> pow(const complex<mpreal>& a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const complex<mpreal>& a, const int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const complex<mpreal>& a, const long double b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const complex<mpreal>& a, const double b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const unsigned int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long double a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const double a, const complex<mpreal>& b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const unsigned long int a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned long int a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned long int a, const int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned long int a, const long double b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned long int a, const double b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const unsigned int a, const unsigned long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned int a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned int a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned int a, const int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned int a, const long double b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const unsigned int a, const double b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const long int a, const unsigned long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const long double b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long int a, const double b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const int a, const unsigned long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const int a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const int a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const int a, const int b, mpc_rnd_t rnd_mode); 
const complex<mpreal> pow(const int a, const long double b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const int a, const double b, mpc_rnd_t rnd_mode); 

const complex<mpreal> pow(const long double a, const long double b, mpc_rnd_t rnd_mode);    
const complex<mpreal> pow(const long double a, const unsigned long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long double a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long double a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const long double a, const int b, mpc_rnd_t rnd_mode);

const complex<mpreal> pow(const double a, const double b, mpc_rnd_t rnd_mode);    
const complex<mpreal> pow(const double a, const unsigned long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const double a, const unsigned int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const double a, const long int b, mpc_rnd_t rnd_mode);
const complex<mpreal> pow(const double a, const int b, mpc_rnd_t rnd_mode);

inline const complex<mpreal> mul_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode);
inline const complex<mpreal> mul_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode);
inline const complex<mpreal> div_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode);
inline const complex<mpreal> div_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode);

//////////////////////////////////////////////////////////////////////////
// 'Dirty' equality check 1: |a-b| < min{|a|,|b|} * eps
inline bool isEqualFuzzy(const complex<mpreal>& a, const complex<mpreal>& b, const mpreal& eps);

// 'Dirty' equality check 2: |a-b| < min{|a|,|b|} * eps( min{|a|,|b|} )
inline bool isEqualFuzzy(const complex<mpreal>& a, const complex<mpreal>& b);

// 'Dirty' realness check: |im(x)| < |x| * eps(|x|)
inline bool isRealFuzzy(const complex<mpreal> x);

// 'Bitwise' equality check
//  maxUlps - a and b can be apart by maxUlps binary numbers. 
inline bool isEqualUlps(const complex<mpreal>& a, const complex<mpreal>& b, int maxUlps);

}
//////////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////////

namespace std {*/
//////////////////////////////////////////////////////////////////////////
// Operators - Assignment
inline complex<mpreal>& complex<mpreal>::operator=(const complex<mpreal>& v)
{
    if (this != &v)
    {
        mp_prec_t tp = mpc_get_prec(  mpc_srcptr());
        mp_prec_t vp = mpc_get_prec(v.mpc_srcptr());

        if(tp != vp){
            clear(mpc_ptr());
            mpc_init2(mpc_ptr(), vp);
        }

        mpc_set(mpc_ptr(), v.mpc_srcptr(), complex<mpreal>::get_default_rnd());

        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    return *this;
}

template <typename real_t> 
inline complex<mpreal>& complex<mpreal>::operator= (const std::complex<real_t>& z)
{
    real(z.real());
    imag(z.imag());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<typename T>
inline complex<mpreal>& complex<mpreal>::operator=(const T v)
{
    mpreal u(v);
    mpc_set_fr(mpc_ptr(), u.mpfr_srcptr(), get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<typename T>
inline complex<mpreal>& complex<mpreal>::real(const T v)
{
    mpreal u(v);
    mpfr_set(mpc_re_ptr(), u.mpfr_srcptr(), MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<typename T>
inline complex<mpreal>& complex<mpreal>::imag(const T v)
{
    mpreal u(v);
    mpfr_set(mpc_im_ptr(), u.mpfr_srcptr(), MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const mpreal& v)
{
    mp_prec_t tp = mpc_get_prec (  mpc_srcptr() );
    mp_prec_t vp = mpfr_get_prec(v.mpfr_srcptr());

    if(tp != vp){
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
    }

    mpc_set_fr(mpc_ptr(), v.mpfr_srcptr(), get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const mpreal& v)
{
    mp_prec_t tp = mpc_get_prec (  mpc_srcptr() );
    mp_prec_t vp = mpfr_get_prec(v.mpfr_srcptr());

    if(tp != vp){
        mpreal im(mpc_im_srcptr());
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), v.mpfr_srcptr(), im.mpfr_srcptr(), get_default_rnd());
    } else {
        mpfr_set(mpc_re_ptr(), v.mpfr_srcptr(), MPC_RND_RE(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const mpreal& v)
{
    mp_prec_t tp = mpc_get_prec (  mpc_srcptr() );
    mp_prec_t vp = mpfr_get_prec(v.mpfr_srcptr());

    if(tp != vp){
        mpreal re(mpc_im_srcptr());
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), v.mpfr_srcptr(), get_default_rnd());
    } else {
        mpfr_set(mpc_im_ptr(), v.mpfr_srcptr(), MPC_RND_IM(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const unsigned int v)        
{    
    mpc_set_ui(mpc_ptr(), v, get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const unsigned int v)
{
    mpfr_set_ui(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const unsigned int v)
{
    mpfr_set_ui(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const unsigned long int v)    
{    
    mpc_set_ui(mpc_ptr(), v, get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const unsigned long int v)
{
    mpfr_set_ui(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const unsigned long int v)
{
    mpfr_set_ui(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const int v)
{    
    mpc_set_si(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const int v)
{
    mpfr_set_si(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const int v)
{
    mpfr_set_si(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const long int v)            
{    
    mpc_set_si(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const long int v)
{
    mpfr_set_si(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const long int v)
{
    mpfr_set_si(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const unsigned long long int v)    
{    
    mpc_set_uj(mpc_ptr(), v, get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const unsigned long long int v)
{
    mpfr_set_uj(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const unsigned long long int v)
{
    mpfr_set_uj(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const long long int v)    
{    
    mpc_set_sj(mpc_ptr(), v, get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const long long int v)
{
    mpfr_set_sj(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const long long int v)
{
    mpfr_set_sj(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const double v)                
{   
#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(v, MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpc_set_d(mpc_ptr(),v,get_default_rnd());
    }else
        throw conversion_overflow();
#else
    mpc_set_d(mpc_ptr(),v,get_default_rnd());
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const double v)
{
#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(v, MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpfr_set_d(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));
    }else
        throw conversion_overflow();
#else
    mpfr_set_d(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const double v)
{
#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
    if(mpfr::mpreal::fits_in_bits(v, MPREAL_DOUBLE_BITS_OVERFLOW))
    {
        mpfr_set_d(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));
    }else
        throw conversion_overflow();
#else
    mpfr_set_d(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const long double v)        
{    
    mpc_set_ld(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const long double v)
{
    mpfr_set_ld(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const long double v)
{
    mpfr_set_ld(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const mpz_t& v)
{
    mpc_set_z(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const mpz_t& v)
{
    mpfr_set_z(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const mpz_t& v)
{
    mpfr_set_z(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const mpq_t& v)
{
    mpc_set_q(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const mpq_t& v)
{
    mpfr_set_q(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const mpq_t& v)
{
    mpfr_set_q(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const mpf_t& v)
{
    mp_prec_t tp = mpc_get_prec(mpc_srcptr());
    mp_prec_t vp = mpf_get_prec(v);

    if(tp != vp){
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
    }

    mpc_set_f(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const mpf_t& v)
{
    mp_prec_t tp = mpc_get_prec(mpc_srcptr());
    mp_prec_t vp = mpf_get_prec(v);

    if(tp != vp){
        mpreal re(v);
        mpreal im(mpc_im_srcptr());
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), im.mpfr_srcptr(), get_default_rnd());
    } else {
        mpfr_set_f(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const mpf_t& v)
{
    mp_prec_t tp = mpc_get_prec(mpc_srcptr());
    mp_prec_t vp = mpf_get_prec (v);

    if(tp != vp){
        mpreal re(mpc_im_srcptr());
        mpreal im(v);
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), im.mpfr_srcptr(), get_default_rnd());
    } else {
        mpfr_set_f(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::operator=(const mpfr_t& v)
{
    mp_prec_t tp = mpc_get_prec (mpc_srcptr());
    mp_prec_t vp = mpfr_get_prec(v);

    if(tp != vp){
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
    }

    mpc_set_fr(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::real(const mpfr_t& v)
{
    mp_prec_t tp = mpc_get_prec (mpc_srcptr());
    mp_prec_t vp = mpfr_get_prec(v);

    if(tp != vp){
        mpreal im(mpc_im_srcptr());
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), v, im.mpfr_srcptr(), get_default_rnd());
    } else {
        mpfr_set(mpc_re_ptr(), v, MPC_RND_RE(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal>& complex<mpreal>::imag(const mpfr_t& v)
{
    mp_prec_t tp = mpc_get_prec (mpc_srcptr());
    mp_prec_t vp = mpfr_get_prec(v);

    if(tp != vp){
        mpreal re(mpc_im_srcptr());
        clear(mpc_ptr());
        mpc_init2(mpc_ptr(), vp);
        mpc_set_fr_fr(mpc_ptr(), re.mpfr_srcptr(), v, get_default_rnd());
    } else {
        mpfr_set(mpc_im_ptr(), v, MPC_RND_IM(get_default_rnd()));
    }

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

#ifdef _COMPLEX_H
inline complex<mpreal>& complex<mpreal>::operator=(const double _Complex v){
    mpc_set_dc(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator=(const long double _Complex v){
    mpc_set_ldc(mpc_ptr(), v, get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}
#endif

inline complex<mpreal>& complex<mpreal>::operator=(const char* s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpcomplex(const char* s,        mp_prec_t prec, int base, mpc_rnd_t mode)
    //        mpcomplex(const std::string& s,mp_prec_t prec, int base, mpc_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpc_t t;

    mpc_init2(t, mpc_get_prec(mpc_srcptr()));

    if(0 == mpc_set_str(t, s, 10, get_default_rnd()))
    {
        mpc_set(mpc_ptr(), t, get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}

inline complex<mpreal>& complex<mpreal>::real(const char* s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s, 10, MPC_RND_RE(get_default_rnd())))
    {
        mpfr_set(mpc_re_ptr(), t, MPC_RND_RE(get_default_rnd()));
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    mpfr_clear(t);
    return *this;
}

inline complex<mpreal>& complex<mpreal>::imag(const char* s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s, 10, MPC_RND_IM(get_default_rnd())))
    {
        mpfr_set(mpc_im_ptr(), t, MPC_RND_IM(get_default_rnd()));
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    mpfr_clear(t);
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator=(const std::string& s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpcomplex(const char* s,        mp_prec_t prec, int base, mpc_rnd_t mode)
    //        mpcomplex(const std::string& s,mp_prec_t prec, int base, mpc_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpc_t t;

    mpc_init2(t, mpc_get_prec(mpc_srcptr()));

    if(0 == mpc_set_str(t, s.c_str(), 10, get_default_rnd()))
    {
        mpc_set(mpc_ptr(), t, get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}

inline complex<mpreal>& complex<mpreal>::real(const std::string& s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s.c_str(), 10, MPC_RND_RE(get_default_rnd())))
    {
        mpfr_set(mpc_re_ptr(), t, MPC_RND_RE(get_default_rnd()));
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    mpfr_clear(t);
    return *this;
}

inline complex<mpreal>& complex<mpreal>::imag(const std::string& s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s.c_str(), 10, MPC_RND_IM(get_default_rnd())))
    {
        mpfr_set(mpc_im_ptr(), t, MPC_RND_IM(get_default_rnd()));
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    mpfr_clear(t);
    return *this;
}

//////////////////////////////////////////////////////////////////////////
// + Addition
inline complex<mpreal>& complex<mpreal>::operator+=(const complex<mpreal>& v)
{
    mpc_add(mpc_ptr(), mpc_srcptr(), v.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const mpreal& v)
{
    mpc_add_fr(mpc_ptr(), mpc_srcptr(), v.mpfr_srcptr(), complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const mpf_t u)
{
    *this += complex<mpreal>(u);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+= (const long double u)
{
    *this += complex<mpreal>(u);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline complex<mpreal>& complex<mpreal>::operator+= (const double u)
{
    *this += complex<mpreal>(u);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const unsigned long int u)
{
    mpc_add_ui(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const unsigned int u)
{
    mpc_add_ui(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const long int u)
{
    mpc_add_si(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const int u)
{
    mpc_add_si(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator+=(const long long int u)         {    *this += complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator+=(const unsigned long long int u){    *this += complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator-=(const long long int  u)        {    *this -= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator-=(const unsigned long long int u){    *this -= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator*=(const long long int  u)        {    *this *= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator*=(const unsigned long long int u){    *this *= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator/=(const long long int  u)        {    *this /= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline complex<mpreal>& complex<mpreal>::operator/=(const unsigned long long int u){    *this /= complex<mpreal>(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }

template<>
inline complex<mpreal> operator+(const complex<mpreal>& a, const complex<mpreal>& b)
{
    complex<mpreal> c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
    mpc_add(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return c;
}

template<> inline complex<mpreal> operator+(const complex<mpreal>& a, const mpreal& b) { return complex<mpreal>(a) += b; }
template<> inline complex<mpreal> operator+(const mpreal& a, const complex<mpreal>& b) { return complex<mpreal>(b) += a; }

inline complex<mpreal>& complex<mpreal>::operator++() 
{
    return *this += 1;
}

inline const complex<mpreal> complex<mpreal>::operator++ (int)
{
    complex<mpreal> x(*this);
    *this += 1;
    return x;
}

inline complex<mpreal>& complex<mpreal>::operator--() 
{
    return *this -= 1;
}

inline const complex<mpreal> complex<mpreal>::operator-- (int)
{
    complex<mpreal> x(*this);
    *this -= 1;
    return x;
}

//////////////////////////////////////////////////////////////////////////
// - Subtraction
inline complex<mpreal>& complex<mpreal>::operator-=(const complex<mpreal>& v)
{
    mpc_sub(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const mpreal& v)
{
    mpc_sub_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const long double v)
{
    *this -= complex<mpreal>(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline complex<mpreal>& complex<mpreal>::operator-=(const double v)
{
    *this -= complex<mpreal>(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const unsigned long int v)
{
    mpc_sub_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const unsigned int v)
{
    mpc_sub_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const long int v)
{
    mpreal V(v);
    mpc_sub_fr(mpc_ptr(),mpc_srcptr(),V.mpfr_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator-=(const int v)
{
    mpreal V(v);
    mpc_sub_fr(mpc_ptr(),mpc_srcptr(),V.mpfr_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const complex<mpreal> complex<mpreal>::operator-()const
{
    complex<mpreal> u(*this);
    mpc_neg(u.mpc_ptr(),u.mpc_srcptr(),complex<mpreal>::get_default_rnd());
    return u;
}

template<>
inline complex<mpreal> operator-(const complex<mpreal>& a, const complex<mpreal>& b)
{
    complex<mpreal> c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
    mpc_sub(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return c;
}

template<> inline complex<mpreal> operator-(const complex<mpreal>& a, const mpreal& b) { return complex<mpreal>(a) -= b; }

template<>
inline complex<mpreal> operator-(const mpreal& b, const complex<mpreal>& a)
{
    complex<mpreal> x(a, mpc_get_prec(a.mpc_ptr()));
    mpc_fr_sub(x.mpc_ptr(), b.mpfr_srcptr(), x.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator-(const unsigned long int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_ptr()));
    mpc_ui_sub(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator-(const unsigned int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_ptr()));
    mpc_ui_sub(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator-(const unsigned short int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_ptr()));
    mpc_ui_sub(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

//////////////////////////////////////////////////////////////////////////
// * Multiplication
inline complex<mpreal>& complex<mpreal>::operator*= (const complex<mpreal>& v)
{
    mpc_mul(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const mpreal& v)
{
    mpc_mul_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const long double v)
{
    *this *= complex<mpreal>(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline complex<mpreal>& complex<mpreal>::operator*=(const double v)
{
    *this *= complex<mpreal>(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const unsigned long int v)
{
    mpc_mul_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const unsigned int v)
{
    mpc_mul_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const long int v)
{
    mpc_mul_si(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator*=(const int v)
{
    mpc_mul_si(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal> operator*(const complex<mpreal>& a, const complex<mpreal>& b)
{
    complex<mpreal> c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
    mpc_mul(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return c;
}

template<> inline complex<mpreal> operator*(const complex<mpreal>& a, const mpreal& b) { return complex<mpreal>(a) *= b; }
template<> inline complex<mpreal> operator*(const mpreal& a, const complex<mpreal>& b) { return complex<mpreal>(b) *= a; }

//////////////////////////////////////////////////////////////////////////
// / Division
inline complex<mpreal>& complex<mpreal>::operator/=(const complex<mpreal>& v)
{
    mpc_div(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const mpreal& v)
{
    mpc_div_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const long double v)
{
    *this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline complex<mpreal>& complex<mpreal>::operator/=(const double v)
{
    *this /= mpreal(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const unsigned long int v)
{
    mpc_div_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const unsigned int v)
{
    mpc_div_ui(mpc_ptr(),mpc_srcptr(),v,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const long int v)
{
    *this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator/=(const int v)
{
    *this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

template<>
inline complex<mpreal> operator/(const complex<mpreal>& a, const complex<mpreal>& b)
{
    complex<mpreal> c(a, (std::max)(mpc_get_prec(a.mpc_srcptr()), mpc_get_prec(b.mpc_srcptr())));
    mpc_div(c.mpc_ptr(), c.mpc_srcptr(), b.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return c;
}

template<> inline complex<mpreal> operator/(const complex<mpreal>& a, const mpreal& b) { return complex<mpreal>(a) /= b; }

template<>
inline complex<mpreal> operator/(const mpreal& b, const complex<mpreal>& a)
{
    complex<mpreal> x(a, mpc_get_prec(a.mpc_ptr()));
    mpc_fr_div(x.mpc_ptr(), b.mpfr_srcptr(), x.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator/(const unsigned long int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_srcptr()));
    mpc_ui_div(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator/(const unsigned int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_srcptr()));
    mpc_ui_div(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

template<>
inline complex<mpreal> operator/(const unsigned short int& b, const complex<mpreal>& a)
{
    complex<mpreal> x(0, mpc_get_prec(a.mpc_srcptr()));
    mpc_ui_div(x.mpc_ptr(), b, a.mpc_srcptr(), complex<mpreal>::get_default_rnd());
    return x;
}

//////////////////////////////////////////////////////////////////////////
// Shifts operators - Multiplication/Division by power of 2
inline complex<mpreal>& complex<mpreal>::operator<<=(const unsigned long int u)
{
    mpc_mul_2ui(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator<<=(const unsigned int u)
{
    mpc_mul_2ui(mpc_ptr(),mpc_srcptr(),static_cast<unsigned long int>(u),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator<<=(const long int u)
{
    mpc_mul_2si(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator<<=(const int u)
{
    mpc_mul_2si(mpc_ptr(),mpc_srcptr(),static_cast<long int>(u),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator>>=(const unsigned long int u)
{
    mpc_div_2ui(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator>>=(const unsigned int u)
{
    mpc_div_2ui(mpc_ptr(),mpc_srcptr(),static_cast<unsigned long int>(u),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator>>=(const long int u)
{
    mpc_div_2si(mpc_ptr(),mpc_srcptr(),u,complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::operator>>=(const int u)
{
    mpc_div_2si(mpc_ptr(),mpc_srcptr(),static_cast<long int>(u),complex<mpreal>::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const complex<mpreal> operator<<(const complex<mpreal>& v, const unsigned long int k)
{
    return mpfr::mul_2ui(v,k);
}

inline const complex<mpreal> operator<<(const complex<mpreal>& v, const unsigned int k)
{
    return mpfr::mul_2ui(v,static_cast<unsigned long int>(k));
}

inline const complex<mpreal> operator<<(const complex<mpreal>& v, const long int k)
{
    return mpfr::mul_2si(v,k);
}

inline const complex<mpreal> operator<<(const complex<mpreal>& v, const int k)
{
    return mpfr::mul_2si(v,static_cast<long int>(k));
}

inline const complex<mpreal> operator>>(const complex<mpreal>& v, const unsigned long int k)
{
    return mpfr::div_2ui(v,k);
}

inline const complex<mpreal> operator>>(const complex<mpreal>& v, const long int k)
{
    return mpfr::div_2si(v,k);
}

inline const complex<mpreal> operator>>(const complex<mpreal>& v, const unsigned int k)
{
    return mpfr::div_2ui(v,static_cast<unsigned long int>(k));
}

inline const complex<mpreal> operator>>(const complex<mpreal>& v, const int k)
{
    return mpfr::div_2si(v,static_cast<long int>(k));
}

//////////////////////////////////////////////////////////////////////////
//Relational operators
template<> inline bool operator == (const complex<mpreal>& a, const complex<mpreal>& b        ){  return (mpc_cmp(a.mpc_srcptr(),b.mpc_srcptr()) != 0 );                   }
inline bool operator == (const complex<mpreal>& a, const unsigned long int b ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const complex<mpreal>& a, const unsigned int b      ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const complex<mpreal>& a, const long int b          ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const complex<mpreal>& a, const int b               ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }

inline bool operator != (const complex<mpreal>& a, const complex<mpreal>& b        ){  return !(a == b);  }
inline bool operator != (const complex<mpreal>& a, const unsigned long int b ){  return !(a == b);  }
inline bool operator != (const complex<mpreal>& a, const unsigned int b      ){  return !(a == b);  }
inline bool operator != (const complex<mpreal>& a, const long int b          ){  return !(a == b);  }
inline bool operator != (const complex<mpreal>& a, const int b               ){  return !(a == b);  }

//////////////////////////////////////////////////////////////////////////
// Type Converters
#ifdef _COMPLEX_H
inline _Complex float complex<mpreal>::toCFloat (mpc_rnd_t mode) const
{
    _Complex float x = (_Complex float) mpc_get_dc(mpc_srcptr(), mode);
    return x;
}

inline std::complex<float> complex<mpreal>::toFloat (mpc_rnd_t mode) const
{
    _Complex float n = (_Complex float) mpc_get_dc(mpc_srcptr(), mode);
    std::complex<float> nxx = std::complex<float>(creal(n), cimag(n));
    return nxx;
}

inline _Complex double complex<mpreal>::toCDouble (mpc_rnd_t mode) const
{
    _Complex double x = mpc_get_dc(mpc_srcptr(), mode);
    return x;
}

inline std::complex<double> complex<mpreal>::toDouble (mpc_rnd_t mode) const
{
    _Complex double n = mpc_get_dc(mpc_srcptr(), mode);
    std::complex<double> nxx = std::complex<double>(creal(n), cimag(n));
    return nxx;
}

inline _Complex long double complex<mpreal>::toCLDouble (mpc_rnd_t mode) const
{
    _Complex long double x = mpc_get_ldc(mpc_srcptr(), mode);
    return x;
}

inline std::complex<long double> complex<mpreal>::toLDouble (mpc_rnd_t mode) const
{
    _Complex long double n = mpc_get_ldc(mpc_srcptr(), mode);
    std::complex<long double> nxx = std::complex<long double>(creal(n), cimag(n));
    return nxx;
}
#endif

inline ::mpc_ptr     complex<mpreal>::mpc_ptr()             { return mp; }
inline ::mpc_srcptr  complex<mpreal>::mpc_ptr()    const    { return mp; }
inline ::mpc_srcptr  complex<mpreal>::mpc_srcptr() const    { return mp; }

inline std::string complex<mpreal>::mpfrString(const ::mpfr_srcptr n, const std::string& format) const
{
    char* s = NULL;
    std::string out;
    if(!format.empty())
    {
        if(!(mpfr_asprintf(&s, format.c_str(), n) < 0))
        {
            out.append(s);
            mpfr_free_str(s);
        }
    }
    return out;
}

// this uses the mpfr asprintf for consistency with mpreal
inline std::string complex<mpreal>::toString(int n, int b, mpc_rnd_t mode) const 
{
    if(b == 10){
        std::string out("(");
        std::ostringstream format;
        int digits = (n >= 0) ? n : 1 + mpfr::bits2digits(mpfr_get_prec(mpfr_srcptr()));
        format << "%." << digits << "RNg";
        out += mpfrString(mpc_re_srcptr(), format.str()) + " ";
        out += mpfrString(mpc_im_srcptr(), format.str()) + ")";
        return out;
    } else {
        char* cstr = mpc_get_str(b, n, mpc_srcptr(), mode);
        std::string out(cstr);
        mpc_free_str(cstr);
        return out;
    }
}

/* the version below uses only inbuilt mpc_get_str but it's inconsistent with mpreal
inline std::string mpcomplex::toString(int n, int b, mpc_rnd_t mode) const
{
    char* cstr = mpc_get_str(b, n, mpc_srcptr(), mode);
    std::string out(cstr);
    mpc_free_str(cstr);

    return out;
}*/

//////////////////////////////////////////////////////////////////////////
// I/O
inline std::ostream& complex<mpreal>::output(std::ostream& os) const 
{
    std::ostringstream format;
    const std::ios::fmtflags flags = os.flags();

    format << ((flags & std::ios::showpos) ? "%+" : "%");
    if (os.precision() >= 0)
        format << '.' << os.precision() << "R*"
               << ((flags & std::ios::floatfield) == std::ios::fixed ? 'f' :
                   (flags & std::ios::floatfield) == std::ios::scientific ? 'e' :
                   'g');
    else
        format << "R*e";

    char *s = NULL;
    if(!(mpfr_asprintf(&s, format.str().c_str(),
                        mpfr::mpreal::get_default_rnd(),
                        mpc_re_srcptr())
        < 0))
    {
        os << "(" << std::string(s);
        mpfr_free_str(s);
        if(!(mpfr_asprintf(&s, format.str().c_str(),
                        mpfr::mpreal::get_default_rnd(),
                        mpc_im_srcptr())
        < 0))
        {
            os << " " << std::string(s) << ")";
        }
        mpfr_free_str(s);
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const complex<mpreal>& v)
{
    return v.output(os);
}

inline std::istream& operator>>(std::istream &is, complex<mpreal>& v)
{
    // TODO: use cout::hexfloat and other flags to setup base
    std::string tmp;
    is >> tmp;
    mpc_set_str(v.mpc_ptr(), tmp.c_str(), 10, complex<mpreal>::get_default_rnd());
    return is;
}

//////////////////////////////////////////////////////////////////////////
// Set/Get number properties
inline int complex<mpreal>::getPrecision() const
{
    return int(mpc_get_prec(mpc_srcptr()));
}

inline complex<mpreal>& complex<mpreal>::setPrecision(int Precision, mpc_rnd_t RoundingMode)
{
    mpc_t v;
    mpc_init2(v, Precision);
    mpc_set(v, mpc_srcptr(), RoundingMode);
    mpc_swap(v, mpc_ptr());
    clear(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mp_prec_t complex<mpreal>::get_prec() const
{
    return mpc_get_prec(mpc_srcptr());
}

inline void complex<mpreal>::set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpc_t v;
    mpc_init2(v, prec);
    mpc_set(v, mpc_srcptr(), rnd_mode);
    mpc_swap(v, mpc_ptr());
    clear(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline complex<mpreal>& complex<mpreal>::setInf(int ReSign, int ImSign)
{
    if(ReSign || ImSign){
        if(ReSign) mpfr_set_inf(mpc_re_ptr(), ReSign);
        if(ImSign) mpfr_set_inf(mpc_im_ptr(), ImSign);
    } else {
        mpfr_set_inf(mpc_re_ptr(), 1);
        mpfr_set_zero(mpc_im_ptr(), 1);
    }
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::setNan() 
{
    mpc_set_nan(mpc_ptr());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::setZero(int ReSign, int ImSign)
{
    mpc_set_si(mpc_ptr(), 0, get_default_rnd());
    setSigns(ReSign, ImSign, get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline complex<mpreal>& complex<mpreal>::setSigns(int ReSign, int ImSign, mpc_rnd_t rnd_mode)
{
    mpfr_setsign(mpc_re_ptr(), mpc_re_srcptr(), ReSign, MPC_RND_RE(rnd_mode));
    mpfr_setsign(mpc_im_ptr(), mpc_im_srcptr(), ImSign, MPC_RND_IM(rnd_mode));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpreal complex<mpreal>::real() const { return mpfr::real(*this, get_default_rnd()); }
inline mpreal complex<mpreal>::imag() const { return mpfr::imag(*this, get_default_rnd()); }

} // end of namespace std

//////////////////////////////////////////////////////////////////////////
// mpfr Namespace Mathematical Functions
//////////////////////////////////////////////////////////////////////////

namespace mpfr {
using std::complex;

//////////////////////////////////////////////////////////////////////////
// Functions needed for std::complex specializations
inline const mpreal abs (const complex<mpreal>& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_abs(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal norm (const complex<mpreal>& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_norm(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal real (const complex<mpreal>& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_real(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal imag (const complex<mpreal>& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_imag(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal arg  (const complex<mpreal>& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_arg(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const complex<mpreal> conj (const complex<mpreal>& x, mpc_rnd_t r)
{
    complex<mpreal> y(x);
    mpc_conj(y.mpc_ptr(), x.mpc_srcptr(), r);
    return y;
}

inline const complex<mpreal> proj (const complex<mpreal>& x, mpc_rnd_t rnd_mode)
{
    complex<mpreal> y(x);
    mpc_proj(y.mpc_ptr(), x.mpc_srcptr(), rnd_mode);
    return y;
}

inline const complex<mpreal> polar(const mpreal& r, const mpreal& theta, mpc_rnd_t rnd_mode)
{
    complex<mpreal> z(r*cos(theta, MPC_RND_RE(rnd_mode)), r*sin(theta, MPC_RND_RE(rnd_mode)));
    return z;
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(a);
    mpc_pow(x.mpc_ptr(),x.mpc_srcptr(),b.mpc_srcptr(),rnd_mode);
    return x;
}
inline const complex<mpreal> pow(const mpreal& a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(a);
    mpc_pow(x.mpc_ptr(),x.mpc_srcptr(),b.mpc_srcptr(),rnd_mode);
    return x;
}
inline const complex<mpreal> pow(const complex<mpreal>& a, const mpreal& b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(b);
    mpc_pow(x.mpc_ptr(),a.mpc_srcptr(),x.mpc_srcptr(),rnd_mode);
    return x;
}

#define MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(f)                    \
        complex<mpreal> y(0, mpc_get_prec(x.mpc_srcptr()));          \
        mpc_##f(y.mpc_ptr(), x.mpc_srcptr(), r);           \
        return y; 

inline const complex<mpreal> sqrt  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqrt );    }
inline const complex<mpreal> log   (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(log  );    }
inline const complex<mpreal> log10 (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(log10);    }
inline const complex<mpreal> exp   (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(exp  );    }
inline const complex<mpreal> cos   (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cos  );    }
inline const complex<mpreal> sin   (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sin  );    }
inline const complex<mpreal> tan   (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tan  );    }
inline const complex<mpreal> acos  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acos );    }
inline const complex<mpreal> asin  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asin );    }
inline const complex<mpreal> atan  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atan );    }
inline const complex<mpreal> cosh  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cosh );    }
inline const complex<mpreal> sinh  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sinh );    }
inline const complex<mpreal> tanh  (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tanh );    }
inline const complex<mpreal> acosh (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acosh);    }
inline const complex<mpreal> asinh (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asinh);    }
inline const complex<mpreal> atanh (const complex<mpreal>& x, mpc_rnd_t r) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atanh);    }

//////////////////////////////////////////////////////////////////////////
// Functions not allowed in std::complex

inline const complex<mpreal> mul_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(v);
    mpc_mul_2ui(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const complex<mpreal> mul_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(v);
    mpc_mul_2si(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const complex<mpreal> div_2ui(const complex<mpreal>& v, unsigned long int k, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(v);
    mpc_div_2ui(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const complex<mpreal> div_2si(const complex<mpreal>& v, long int k, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(v);
    mpc_div_2si(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const complex<mpreal> sqr  (const complex<mpreal>& x, mpc_rnd_t r)
{   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqr );    }

inline int sin_cos(complex<mpreal>& s, complex<mpreal>& c, const complex<mpreal>& v, mpc_rnd_t rnd_mode_s, mpc_rnd_t rnd_mode_c)
{
    return mpc_sin_cos(s.mpc_ptr(), c.mpc_ptr(), v.mpc_srcptr(), rnd_mode_s, rnd_mode_c);
}

inline const complex<mpreal> sqrtcomp (const mpreal v, mpc_rnd_t rnd_mode)               {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const unsigned long int v, mpc_rnd_t rnd_mode)    {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const unsigned int v, mpc_rnd_t rnd_mode)         {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const long int v, mpc_rnd_t rnd_mode)             {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const int v, mpc_rnd_t rnd_mode)                  {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const long double v, mpc_rnd_t rnd_mode)          {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }
inline const complex<mpreal> sqrt     (const double v, mpc_rnd_t rnd_mode)               {   return mpfr::sqrt(complex<mpreal>(v),rnd_mode);    }

inline const complex<mpreal> sec   (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::cos(x, r);                           }
inline const complex<mpreal> csc   (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::sin(x, r);                           }
inline const complex<mpreal> cot   (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::tan(x, r);                           }
inline const complex<mpreal> sech  (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::cosh(x, r);                          }
inline const complex<mpreal> csch  (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::sinh(x, r);                          }
inline const complex<mpreal> coth  (const complex<mpreal>& x, mpc_rnd_t r) {   return 1/mpfr::tanh(x, r);                          }
inline const complex<mpreal> acot  (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::atan (1/v, r);                      }
inline const complex<mpreal> asec  (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::acos (1/v, r);                      }
inline const complex<mpreal> acsc  (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::asin (1/v, r);                      }
inline const complex<mpreal> acoth (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::atanh(1/v, r);                      }
inline const complex<mpreal> asech (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::acosh(1/v, r);                      }
inline const complex<mpreal> acsch (const complex<mpreal>& v, mpc_rnd_t r) {   return mpfr::asinh(1/v, r);                      }

inline const complex<mpreal> fma (const complex<mpreal>& v1, const complex<mpreal>& v2, const complex<mpreal>& v3, mpc_rnd_t rnd_mode)
{
    complex<mpreal> a;
    mp_prec_t p1, p2, p3;

    p1 = v1.get_prec(); 
    p2 = v2.get_prec(); 
    p3 = v3.get_prec(); 

    a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

    mpc_fma(a.mpc_ptr(),v1.mpc_srcptr(),v2.mpc_srcptr(),v3.mpc_srcptr(),rnd_mode);
    return a;
}

//////////////////////////////////////////////////////////////////////////
// mpreal member functions not allowed in namespace std

//////////////////////////////////////////////////////////////////////////
// Miscellaneous Functions

inline void swap(complex<mpreal>& x, complex<mpreal>& y) { mpc_swap(x.mpc_ptr(), y.mpc_ptr()); }

inline const complex<mpreal> mpc_urandom (gmp_randstate_t& state)
{
    complex<mpreal> x;
    mpc_urandom(x.mpc_ptr(), state);
    return x;
}

inline const complex<mpreal> mpc_const_infinity (int ReSign = 1, int ImSign = 0, mp_prec_t p = mpreal::get_default_prec())
{
    complex<mpreal> x(0, p);
    if(ReSign || ImSign){
        if(ReSign) mpfr_set_inf(x.mpc_re_ptr(), ReSign);
        if(ImSign) mpfr_set_inf(x.mpc_im_ptr(), ImSign);
    } else {
        mpfr_set_inf(x.mpc_re_ptr(), 1);
        mpfr_set_zero(x.mpc_im_ptr(), 1);
    }
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return x;
}

inline bool isEqualUlps(const complex<mpreal>& a, const complex<mpreal>& b, int maxUlps)
{
    return mpfr::abs(a - b) <= machine_epsilon((max)(abs(a), abs(b))) * maxUlps;
}

inline bool isEqualFuzzy(const complex<mpreal>& a, const complex<mpreal>& b, const mpreal& eps)
{
    return mpfr::abs(a - b) <= eps;
}

inline bool isEqualFuzzy(const complex<mpreal>& a, const complex<mpreal>& b)
{
    return isEqualFuzzy(a, b, machine_epsilon((max)(mpreal(1), (min)(abs(a), abs(b)))));
}

inline bool isRealFuzzy(const complex<mpreal> x)
{
    return isEqualFuzzy(imag(x), mpreal(0), machine_epsilon(abs(x)));
}

//////////////////////////////////////////////////////////////////////////
// The Great Pow Sea

inline const complex<mpreal> pow(const complex<mpreal>& a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(a);
    mpc_pow_ui(x.mpc_ptr(),x.mpc_srcptr(),b,rnd_mode);
    return x;
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(a,static_cast<unsigned long int>(b),rnd_mode);
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const long int b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(a);
    mpc_pow_si(x.mpc_ptr(),x.mpc_srcptr(),b,rnd_mode);
    return x;
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,static_cast<long int>(b),rnd_mode);
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const long double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,complex<mpreal>(b),rnd_mode);
}

inline const complex<mpreal> pow(const complex<mpreal>& a, const double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,complex<mpreal>(b),rnd_mode);
}

inline const complex<mpreal> pow(const unsigned long int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    complex<mpreal> x(a);
    mpc_pow(x.mpc_ptr(),x.mpc_srcptr(),b.mpc_srcptr(),rnd_mode);
    return x;
}

inline const complex<mpreal> pow(const unsigned int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode);
}

inline const complex<mpreal> pow(const long int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    if(a>=0)  return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode);
    else      return mpfr::pow(complex<mpreal>(a),b,rnd_mode);
}

inline const complex<mpreal> pow(const int a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    if(a>=0)  return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode);
    else      return mpfr::pow(complex<mpreal>(a),b,rnd_mode);
}

inline const complex<mpreal> pow(const long double a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode);
}

inline const complex<mpreal> pow(const double a, const complex<mpreal>& b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode);
}

// pow unsigned long int
inline const complex<mpreal> pow(const unsigned long int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    mpreal x(a);
    mpfr_ui_pow_ui(x.mpfr_ptr(),a,b,MPC_RND_RE(rnd_mode));
    return complex<mpreal>(x);
}

inline const complex<mpreal> pow(const unsigned long int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const complex<mpreal> pow(const unsigned long int a, const long int b, mpc_rnd_t rnd_mode)
{
    if(b>=0)    return mpfr::pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else        return mpfr::pow(a,complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned long int a, const int b, mpc_rnd_t rnd_mode)
{
    if(b>=0)    return mpfr::pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else        return mpfr::pow(a,complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned long int a, const long double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned long int a, const double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(a,complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

// pow unsigned int
inline const complex<mpreal> pow(const unsigned int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
}

inline const complex<mpreal> pow(const unsigned int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const complex<mpreal> pow(const unsigned int a, const long int b, mpc_rnd_t rnd_mode)
{
    if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned int a, const int b, mpc_rnd_t rnd_mode)
{
    if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned int a, const long double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

inline const complex<mpreal> pow(const unsigned int a, const double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
}

// pow long int
inline const complex<mpreal> pow(const long int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    if (a>=0) return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
    else      return mpfr::pow(complex<mpreal>(a),b,rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const long int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    if (a>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
    else      return mpfr::pow(complex<mpreal>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const long int a, const long int b, mpc_rnd_t rnd_mode)
{
    if (a>=0)
    {
        if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    }else{
        return mpfr::pow(complex<mpreal>(a),b,rnd_mode); // mpfr_pow_si
    }
}

inline const complex<mpreal> pow(const long int a, const int b, mpc_rnd_t rnd_mode)
{
    if (a>=0)
    {
        if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    }else{
        return mpfr::pow(complex<mpreal>(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
    }
}

inline const complex<mpreal> pow(const long int a, const long double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    else        return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode); //mpfr_pow
}

inline const complex<mpreal> pow(const long int a, const double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    else        return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode); //mpfr_pow
}

// pow int
inline const complex<mpreal> pow(const int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    if (a>=0) return mpfr::pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
    else      return mpfr::pow(complex<mpreal>(a),b,rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    if (a>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
    else      return mpfr::pow(complex<mpreal>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const int a, const long int b, mpc_rnd_t rnd_mode)
{
    if (a>=0)
    {
        if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    }else{
        return mpfr::pow(complex<mpreal>(a),b,rnd_mode); // mpfr_pow_si
    }
}

inline const complex<mpreal> pow(const int a, const int b, mpc_rnd_t rnd_mode)
{
    if (a>=0)
    {
        if(b>=0) return mpfr::pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else     return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    }else{
        return mpfr::pow(complex<mpreal>(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
    }
}

inline const complex<mpreal> pow(const int a, const long double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    else        return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode); //mpfr_pow
}

inline const complex<mpreal> pow(const int a, const double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return mpfr::pow(static_cast<unsigned long int>(a),complex<mpreal>(b),rnd_mode); //mpfr_ui_pow
    else        return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode); //mpfr_pow
}

// pow long double 
inline const complex<mpreal> pow(const long double a, const long double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode);
}

inline const complex<mpreal> pow(const long double a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const long double a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const complex<mpreal> pow(const long double a, const long int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode); // mpfr_pow_si
}

inline const complex<mpreal> pow(const long double a, const int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}

inline const complex<mpreal> pow(const double a, const double b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),complex<mpreal>(b),rnd_mode);
}

inline const complex<mpreal> pow(const double a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode); // mpfr_pow_ui
}

inline const complex<mpreal> pow(const double a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),static_cast<unsigned long int>(b),rnd_mode); // mpfr_pow_ui
}

inline const complex<mpreal> pow(const double a, const long int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),b,rnd_mode); // mpfr_pow_si
}

inline const complex<mpreal> pow(const double a, const int b, mpc_rnd_t rnd_mode)
{
    return mpfr::pow(complex<mpreal>(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}

} // end of namespace mpfr

// we are allowed to extend namespace std with specializations only
namespace std
{
	// Explicit specialization of std::swap for mpcomplex numbers
	// Thus standard algorithms will use efficient version of swap (due to Koenig lookup)
	// Non-throwing swap C++ idiom: http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-throwing_swap
    template<> inline void swap(complex<mpreal>& x, complex<mpreal>& y) { mpfr::swap(x, y); }

	// std::complex function specializations
	template<> inline complex<mpreal> polar(const mpreal& r, const mpreal& theta) { return mpfr::polar(r, theta, complex<mpreal>::get_default_rnd()); }

	template<> inline complex<mpreal> pow(const complex<mpreal>& a, const complex<mpreal>& b) { return mpfr::pow(a, b, complex<mpreal>::get_default_rnd()); }
	/* OS X clang seems to be missing these two overloads despite 26.4.8
	template<> inline complex<mpreal> pow(const complex<mpreal>& a, const          mpreal& b) { return mpfr::pow(a, b, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> pow(const          mpreal& a, const complex<mpreal>& b) { return mpfr::pow(a, b, complex<mpreal>::get_default_rnd()); } */

	template<> inline complex<mpreal> sqrt  (const complex<mpreal>& z) { return mpfr::sqrt(z, complex<mpreal>::get_default_rnd()); }
	template<> inline          mpreal abs   (const complex<mpreal>& z) { return mpfr::abs  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline          mpreal norm  (const complex<mpreal>& z) { return mpfr::norm (z, complex<mpreal>::get_default_rnd()); }
	template<> inline          mpreal real  (const complex<mpreal>& z) { return mpfr::real (z, complex<mpreal>::get_default_rnd()); }
	template<> inline          mpreal imag  (const complex<mpreal>& z) { return mpfr::imag (z, complex<mpreal>::get_default_rnd()); }
	template<> inline          mpreal arg   (const complex<mpreal>& z) { return mpfr::arg  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> conj  (const complex<mpreal>& z) { return mpfr::conj (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> proj  (const complex<mpreal>& z) { return mpfr::proj (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> log   (const complex<mpreal>& z) { return mpfr::log  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> log10 (const complex<mpreal>& z) { return mpfr::log10(z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> exp   (const complex<mpreal>& z) { return mpfr::exp  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> cos   (const complex<mpreal>& z) { return mpfr::cos  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> sin   (const complex<mpreal>& z) { return mpfr::sin  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> tan   (const complex<mpreal>& z) { return mpfr::tan  (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> acos  (const complex<mpreal>& z) { return mpfr::acos (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> asin  (const complex<mpreal>& z) { return mpfr::asin (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> atan  (const complex<mpreal>& z) { return mpfr::atan (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> cosh  (const complex<mpreal>& z) { return mpfr::cosh (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> sinh  (const complex<mpreal>& z) { return mpfr::sinh (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> tanh  (const complex<mpreal>& z) { return mpfr::tanh (z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> acosh (const complex<mpreal>& z) { return mpfr::acosh(z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> asinh (const complex<mpreal>& z) { return mpfr::asinh(z, complex<mpreal>::get_default_rnd()); }
	template<> inline complex<mpreal> atanh (const complex<mpreal>& z) { return mpfr::atanh(z, complex<mpreal>::get_default_rnd()); }

	// We are not allowed to specialize std::numeric_limits for std::complex<T>, cf 18.3.2.1/4
}

#endif /* __MPCOMPLEX_H__ */
