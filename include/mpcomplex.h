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

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <cstdint>
#include <complex.h>
#include <complex>
#include <algorithm>

#include "mpreal.h"
#include <mpc.h>

#ifdef MPREAL_HAVE_MOVE_SUPPORT
    #define mpc_is_initialized(x)      (0 != (x)->re->_mpfr_d || 0 != (x)->im->_mpfr_d)
    #define mpc_set_uninitialized(x)   ((x)->re->_mpfr_d = 0; (x)->im->_mpfr_d = 0 )
#endif

namespace mpfr {

class mpcomplex {
private:
    mpc_t mp;
    
public:
    
    // MPC doesn't have default rnd/prec, so we use mpreal's defaults for both components
    inline static mpc_rnd_t   get_default_rnd()   {    return (mpc_rnd_t)MPC_RND(mpreal::get_default_rnd(), mpreal::get_default_rnd());       }
    inline static mpfr_prec_t  get_default_prec() {    return mpreal::get_default_prec(); }

    // Constructors && type conversions
    mpcomplex();
	mpcomplex(const mpreal& rp);
    mpcomplex(const mpcomplex& u);
	mpcomplex(const mpcomplex& u, mpc_rnd_t mode);
    mpcomplex(const mpf_t u);    
    mpcomplex(const double u,                 mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const long double u,            mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const unsigned long long int u, mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const long long int u,          mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const unsigned long int u,      mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const unsigned int u,           mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const long int u,               mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const int u,                    mp_prec_t prec = mpcomplex::get_default_prec(), mpc_rnd_t mode = mpcomplex::get_default_rnd());
    
    // Construct mpcomplex from mpc_t structure.
    // shared = true allows to avoid deep copy, so that mpcomplex and 'u' share the same data & pointers.    
    mpcomplex(const mpc_t  u, bool shared = false);   

    mpcomplex(const char* s,             mp_prec_t prec = mpcomplex::get_default_prec(), int base = 10, mpc_rnd_t mode = mpcomplex::get_default_rnd());
    mpcomplex(const std::string& s,      mp_prec_t prec = mpcomplex::get_default_prec(), int base = 10, mpc_rnd_t mode = mpcomplex::get_default_rnd());

    ~mpcomplex();                           

#ifdef MPREAL_HAVE_MOVE_SUPPORT
    mpcomplex& operator=(mpcomplex&& v);
    mpcomplex(mpcomplex&& u);
#endif

    // Operations
    // =
    // +, -, *, /, ++, --, <<, >> 
    // *=, +=, -=, /=,
    // <, >, ==, <=, >=

    // =
    mpcomplex& operator=(const mpcomplex& v);
    template <typename real_t> mpcomplex& operator= (const std::complex<real_t>& z);
	mpcomplex& operator=(const mpreal& v);
    mpcomplex& operator=(const mpf_t v);
    mpcomplex& operator=(const long double v);
    mpcomplex& operator=(const double v);        
    mpcomplex& operator=(const unsigned long int v);
    mpcomplex& operator=(const unsigned long long int v);
    mpcomplex& operator=(const long long int v);
    mpcomplex& operator=(const unsigned int v);
    mpcomplex& operator=(const long int v);
    mpcomplex& operator=(const int v);
    mpcomplex& operator=(const char* s);
    mpcomplex& operator=(const std::string& s);
	mpcomplex& set_real(const mpreal& v);
	mpcomplex& set_imag(const mpreal& v);
    mpcomplex& set_real(const mpf_t v);
    mpcomplex& set_imag(const mpf_t v);
    mpcomplex& set_real(const long double v);
    mpcomplex& set_imag(const long double v);
    mpcomplex& set_real(const double v);        
    mpcomplex& set_imag(const double v);        
    mpcomplex& set_real(const unsigned long int v);
    mpcomplex& set_imag(const unsigned long int v);
    mpcomplex& set_real(const unsigned long long int v);
    mpcomplex& set_imag(const unsigned long long int v);
    mpcomplex& set_real(const long long int v);
    mpcomplex& set_imag(const long long int v);
	mpcomplex& set_real(const unsigned int v);
	mpcomplex& set_imag(const unsigned int v);
    mpcomplex& set_real(const long int v);
    mpcomplex& set_imag(const long int v);
    mpcomplex& set_real(const int v);
    mpcomplex& set_imag(const int v);
    mpcomplex& set_real(const char* s);
    mpcomplex& set_imag(const char* s);
    mpcomplex& set_real(const std::string& s);
    mpcomplex& set_imag(const std::string& s);

    // +; iplus(u) represents += i*u
    mpcomplex& operator+=(const mpcomplex& v);
	mpcomplex& operator+=(const mpreal &v);
    mpcomplex& operator+=(const mpf_t v);
    mpcomplex& operator+=(const long double u);
    mpcomplex& operator+=(const double u);
    mpcomplex& operator+=(const unsigned long int u);
    mpcomplex& operator+=(const unsigned int u);
    mpcomplex& operator+=(const long int u);
    mpcomplex& operator+=(const int u);
	mpcomplex& iplus(const mpreal &v);
    mpcomplex& iplus(const mpf_t v);
    mpcomplex& iplus(const long double u);
    mpcomplex& iplus(const double u);
    mpcomplex& iplus(const unsigned long int u);
    mpcomplex& iplus(const unsigned int u);
    mpcomplex& iplus(const long int u);
    mpcomplex& iplus(const int u);

    mpcomplex& operator+=(const long long int  u);
    mpcomplex& operator+=(const unsigned long long int u);
    mpcomplex& iplus(const long long int  u);
    mpcomplex& iplus(const unsigned long long int u);
    mpcomplex& operator-=(const long long int  u);
    mpcomplex& operator-=(const unsigned long long int u);
    mpcomplex& iminus(const long long int  u);
    mpcomplex& iminus(const unsigned long long int u);
    mpcomplex& operator*=(const long long int  u);
    mpcomplex& operator*=(const unsigned long long int u);
    mpcomplex& itimes(const long long int  u);
    mpcomplex& itimes(const unsigned long long int u);
    mpcomplex& operator/=(const long long int  u);
    mpcomplex& operator/=(const unsigned long long int u);
    mpcomplex& idiv(const long long int  u);
    mpcomplex& idiv(const unsigned long long int u);

    const mpcomplex operator+() const;
    mpcomplex& operator++ ();
    const mpcomplex  operator++ (int); 

    // -; iminus(u) represents -= i*u
    mpcomplex& operator-=(const mpcomplex& v);
    mpcomplex& operator-=(const mpreal& v);
    mpcomplex& operator-=(const long double u);
    mpcomplex& operator-=(const double u);
    mpcomplex& operator-=(const unsigned long int u);
    mpcomplex& operator-=(const unsigned int u);
    mpcomplex& operator-=(const long int u);
    mpcomplex& operator-=(const int u);
	mpcomplex& iminus(const mpreal &v);
    mpcomplex& iminus(const mpf_t v);
    mpcomplex& iminus(const long double u);
    mpcomplex& iminus(const double u);
    mpcomplex& iminus(const unsigned long int u);
    mpcomplex& iminus(const unsigned int u);
    mpcomplex& iminus(const long int u);
    mpcomplex& iminus(const int u);
    const mpcomplex operator-() const;
    friend const mpcomplex operator-(const unsigned long int b, const mpcomplex& a);
    friend const mpcomplex operator-(const unsigned int b,      const mpcomplex& a);
    friend const mpcomplex operator-(const long int b,          const mpcomplex& a);
    friend const mpcomplex operator-(const int b,               const mpcomplex& a);
    friend const mpcomplex operator-(const double b,            const mpcomplex& a);
    mpcomplex& operator-- ();    
    const mpcomplex  operator-- (int);

    // *; itimes(u) represents *= i*u
    mpcomplex& operator*=(const mpcomplex& v);
	mpcomplex& operator*=(const mpreal& v);
    mpcomplex& operator*=(const long double v);
    mpcomplex& operator*=(const double v);
    mpcomplex& operator*=(const unsigned long int v);
    mpcomplex& operator*=(const unsigned int v);
    mpcomplex& operator*=(const long int v);
    mpcomplex& operator*=(const int v);
	mpcomplex& itimes(const mpreal &v);
    mpcomplex& itimes(const mpf_t v);
    mpcomplex& itimes(const long double u);
    mpcomplex& itimes(const double u);
    mpcomplex& itimes(const unsigned long int u);
    mpcomplex& itimes(const unsigned int u);
    mpcomplex& itimes(const long int u);
    mpcomplex& itimes(const int u);
    
    // /; idiv(u) represents /= i*u
    mpcomplex& operator/=(const mpcomplex& v);
    mpcomplex& operator/=(const mpreal& v);
    mpcomplex& operator/=(const long double v);
    mpcomplex& operator/=(const double v);
    mpcomplex& operator/=(const unsigned long int v);
    mpcomplex& operator/=(const unsigned int v);
    mpcomplex& operator/=(const long int v);
    mpcomplex& operator/=(const int v);
	mpcomplex& idiv(const mpreal &v);
    mpcomplex& idiv(const mpf_t v);
    mpcomplex& idiv(const long double u);
    mpcomplex& idiv(const double u);
    mpcomplex& idiv(const unsigned long int u);
    mpcomplex& idiv(const unsigned int u);
    mpcomplex& idiv(const long int u);
    mpcomplex& idiv(const int u);
    friend const mpcomplex operator/(const unsigned long int b, const mpcomplex& a);
    friend const mpcomplex operator/(const unsigned int b,      const mpcomplex& a);
    friend const mpcomplex operator/(const long int b,          const mpcomplex& a);
    friend const mpcomplex operator/(const int b,               const mpcomplex& a);
    friend const mpcomplex operator/(const double b,            const mpcomplex& a);

    //<<= Fast Multiplication by 2^u
    mpcomplex& operator<<=(const unsigned long int u);
    mpcomplex& operator<<=(const unsigned int u);
    mpcomplex& operator<<=(const long int u);
    mpcomplex& operator<<=(const int u);

    //>>= Fast Division by 2^u
    mpcomplex& operator>>=(const unsigned long int u);
    mpcomplex& operator>>=(const unsigned int u);
    mpcomplex& operator>>=(const long int u);
    mpcomplex& operator>>=(const int u);

    // Type Conversion operators
	_Complex float              toCFloat    (mpc_rnd_t mode = MPC_RNDNN)    const;
	std::complex<float>         toFloat     (mpc_rnd_t mode = MPC_RNDNN)    const;
	_Complex double             toCDouble   (mpc_rnd_t mode = MPC_RNDNN)    const;
	std::complex<double>        toDouble    (mpc_rnd_t mode = MPC_RNDNN)    const;
	_Complex long double        toCLDouble  (mpc_rnd_t mode = MPC_RNDNN)    const;
	std::complex<long double>   toLDouble   (mpc_rnd_t mode = MPC_RNDNN)    const;

#if defined (mpcomplex_HAVE_EXPLICIT_CONVERTERS)
    explicit operator _Complex float             () const { return toCFloat();               }
	explicit operator std::complex<float>        () const { return toFloat();                }
    explicit operator _Complex double            () const { return toCDouble();              }
	explicit operator std::complex<double>       () const { return toDouble();               }
    explicit operator _Complex long double       () const { return toCLDouble();             }
	explicit operator std::complex<long double>  () const { return toLDouble();              }
#endif

    // Get raw pointers so that mpcomplex can be directly used in raw mpc_* functions
    ::mpc_ptr     mpc_ptr();
    ::mpc_srcptr  mpc_ptr()       const;
    ::mpc_srcptr  mpc_srcptr()    const;
	::mpfr_ptr    mpc_re_ptr()			{ return mpc_realref(mpc_ptr());    }
	::mpfr_srcptr mpc_re_ptr()    const { return mpc_realref(mpc_srcptr()); }
	::mpfr_srcptr mpc_re_srcptr() const { return mpc_realref(mpc_srcptr()); }
	::mpfr_ptr    mpc_im_ptr()			{ return mpc_imagref(mpc_ptr());    }
	::mpfr_srcptr mpc_im_ptr()    const { return mpc_imagref(mpc_srcptr()); }
	::mpfr_srcptr mpc_im_srcptr() const { return mpc_imagref(mpc_srcptr()); }

    // Convert mpcomplex to string with n significant digits in base b
    // n = -1 -> convert with the maximum available digits
    std::string toString(int n = -1, int b = 10, mpc_rnd_t mode = mpcomplex::get_default_rnd()) const;

#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
    std::string mpfrString(const ::mpfr_srcptr n, const std::string& format) const;
#endif

    std::ostream& output(std::ostream& os) const;

    // Math Functions
    friend const mpreal  abs (const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal norm (const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal real (const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());
    friend const mpreal imag (const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());
	friend const mpreal arg (const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());

	friend const mpcomplex conj(const mpcomplex& v, mpc_rnd_t rnd_mode = get_default_rnd());

    friend const mpcomplex sqr (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sqrt(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sqrt(const mpreal& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sqrt(const unsigned long int v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex sqrt(const unsigned int v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex sqrt(const long int v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex sqrt(const int v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex sqrt(const long double v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex sqrt(const double v, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
    friend const mpcomplex root(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode);
    friend const mpcomplex pow (const mpcomplex& a, const mpcomplex& b, mpc_rnd_t rnd_mode);
    friend const mpcomplex pow (const mpcomplex& a, const unsigned long int b, mpc_rnd_t rnd_mode);
    friend const mpcomplex pow (const mpcomplex& a, const long int b, mpc_rnd_t rnd_mode);
    friend const mpcomplex pow (const unsigned long int a, const mpcomplex& b, mpc_rnd_t rnd_mode);
    friend const mpcomplex pow (const unsigned long int a, const unsigned long int b, mpc_rnd_t rnd_mode);

    friend inline const mpcomplex mul_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode);
    friend inline const mpcomplex mul_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode);
    friend inline const mpcomplex div_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode);
    friend inline const mpcomplex div_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode);
    
    friend const mpcomplex log  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex log10(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex exp  (const mpcomplex& v, mpc_rnd_t rnd_mode); 

    friend const mpcomplex cos(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sin(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex tan(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sec(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex csc(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex cot(const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend int sin_cos(mpcomplex& s, mpcomplex& c, const mpcomplex& v, mpc_rnd_t rnd_mode);

    friend const mpcomplex acos  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex asin  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex atan  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex acot  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex asec  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex acsc  (const mpcomplex& v, mpc_rnd_t rnd_mode);

    friend const mpcomplex cosh  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sinh  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex tanh  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex sech  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex csch  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex coth  (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex acosh (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex asinh (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex atanh (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex acoth (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex asech (const mpcomplex& v, mpc_rnd_t rnd_mode);
    friend const mpcomplex acsch (const mpcomplex& v, mpc_rnd_t rnd_mode);

    friend const mpcomplex mpc_urandom (gmp_randstate_t& state);     // use gmp_randinit_default() to init state, gmp_randclear() to clear

	// constant infinity value in case someone cares
    friend const mpreal const_infinity(int ReSign, int ImSign, mp_prec_t prec);

    // Output/ Input
    friend std::ostream& operator<<(std::ostream& os, const mpcomplex& v);
    friend std::istream& operator>>(std::istream& is, mpcomplex& v);

    // Set/Get instance properties
    inline mp_prec_t    get_prec() const;
    inline void         set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode = get_default_rnd());    // Change precision with rounding mode
	inline int          getPrecision() const;
	inline mpcomplex&   setPrecision(int Precision, mpc_rnd_t RoundingMode = get_default_rnd());

    // Set mpcomplex to +/- inf, NaN, +/-0
    mpcomplex&        setInf  (int ReSign = +1, int ImSign = 0);    
    mpcomplex&        setNan  ();
    mpcomplex&        setZero (int ReSign = +1, int ImSign = +1);
    mpcomplex&        setSigns(int ReSign, int ImSign, mpc_rnd_t rnd_mode);

    // Inexact conversion from float
    inline bool fits_in_bits(double x, int n);

    // Set/Get global properties
    static void            set_default_prec(mp_prec_t prec);
    static void            set_default_rnd(mpc_rnd_t rnd_mode);

    static mp_exp_t  get_emin (void);
    static mp_exp_t  get_emax (void);
    static mp_exp_t  get_emin_min (void);
    static mp_exp_t  get_emin_max (void);
    static mp_exp_t  get_emax_min (void);
    static mp_exp_t  get_emax_max (void);
    static int       set_emin (mp_exp_t exp);
    static int       set_emax (mp_exp_t exp);

    // Efficient swapping of two mpcomplex values - needed for std algorithms
    friend void swap(mpcomplex& x, mpcomplex& y);
    
private:
    // Human friendly Debug Preview in Visual Studio.
    // Put one of these lines:
    //
    // mpfr::mpcomplex=<DebugView>                              ; Show value only
    // mpfr::mpcomplex=<DebugView>, <mp[0]._mpfr_prec,u>bits    ; Show value & precision
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
inline mpcomplex::mpcomplex() 
{ 
    mpc_init2(mpc_ptr(), mpcomplex::get_default_prec()); 
    mpc_set_ui(mpc_ptr(), 0, mpcomplex::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const mpcomplex& u) 
{
    mpc_init2(mpc_ptr(), mpc_get_prec(u.mpc_srcptr()));
    mpc_set  (mpc_ptr(), u.mpc_srcptr(), mpcomplex::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const mpcomplex& u, mpc_rnd_t mode)
{
    mpc_init2(mpc_ptr(), mpc_get_prec(u.mpc_srcptr()));
	mpc_set  (mpc_ptr(), u.mpc_srcptr(), mode);

	MPREAL_MSVC_DEBUGVIEW_CODE;
}

#ifdef MPREAL_HAVE_MOVE_SUPPORT
inline mpcomplex::mpcomplex(mpcomplex&& other)
{
    mpfr_set_uninitialized(mpfr_ptr());     // make sure "other" holds no pointer to actual data 
    mpc_swap(mpc_ptr(), other.mpc_ptr());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex& mpcomplex::operator=(mpcomplex&& other)
{
    mpc_swap(mpc_ptr(), other.mpc_ptr());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}
#endif

inline mpcomplex::mpcomplex(const mpc_t  u, bool shared)
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

inline mpcomplex::mpcomplex(const mpreal& rp)
{
    mpc_init2(mpc_ptr(), mpfr_get_prec(rp.mpfr_srcptr()));
    mpc_set_fr(mpc_ptr(), rp.mpfr_srcptr(), get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const mpf_t u)
{
    mpc_init2(mpc_ptr(),(mp_prec_t) mpf_get_prec(u)); // (gmp: mp_bitcnt_t) unsigned long -> long (mpfr: mp_prec_t)
    mpc_set_f(mpc_ptr(),u,mpcomplex::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const double u, mp_prec_t prec, mpc_rnd_t mode)
{
     mpc_init2(mpc_ptr(), prec);

#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
	if(fits_in_bits(u, MPREAL_DOUBLE_BITS_OVERFLOW))
	{
		mpc_set_d(mpc_ptr(), u, mode);
	}else
		throw conversion_overflow();
#else
	mpc_set_d(mpc_ptr(), u, mode);
#endif

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const long double u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_ld(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const unsigned long long int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_uj(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const long long int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_sj(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const unsigned long int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_ui(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const unsigned int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_ui(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const long int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_si(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const int u, mp_prec_t prec, mpc_rnd_t mode)
{ 
    mpc_init2 (mpc_ptr(), prec);
    mpc_set_si(mpc_ptr(), u, mode);

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const char* s, mp_prec_t prec, int base, mpc_rnd_t mode)
{
    mpc_init2  (mpc_ptr(), prec);
    mpc_set_str(mpc_ptr(), s, base, mode); 

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex::mpcomplex(const std::string& s, mp_prec_t prec, int base, mpc_rnd_t mode)
{
    mpc_init2  (mpc_ptr(), prec);
    mpc_set_str(mpc_ptr(), s.c_str(), base, mode); 

    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline void mpcomplex::clear(::mpc_ptr x)
{
#ifdef MPREAL_HAVE_MOVE_SUPPORT
    if(mpc_is_initialized(x)) 
#endif
    mpc_clear(x);
}

inline mpcomplex::~mpcomplex() 
{ 
    clear(mpc_ptr());
}                           

// internal namespace needed for template magic
namespace infernal{

    // Use SFINAE to restrict arithmetic operations instantiation only for numeric types
    // This is needed for smooth integration with libraries based on expression templates, like Eigen.
    // TODO: Do the same for boolean operators.
    template <typename ArgumentType> struct result_type {};    
    
    template <> struct result_type<mpcomplex>           {typedef mpcomplex type;};    
    template <> struct result_type<long double>         {typedef mpcomplex type;};    
    template <> struct result_type<double>              {typedef mpcomplex type;};    
    template <> struct result_type<unsigned long int>   {typedef mpcomplex type;};    
    template <> struct result_type<unsigned int>        {typedef mpcomplex type;};    
    template <> struct result_type<long int>            {typedef mpcomplex type;};    
    template <> struct result_type<int>                 {typedef mpcomplex type;};    
    template <> struct result_type<long long>           {typedef mpcomplex type;};    
    template <> struct result_type<unsigned long long>  {typedef mpcomplex type;};    
}

// + Addition
template <typename Rhs> 
inline const typename infernal::result_type<Rhs>::type 
    operator+(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) += rhs;    }

template <typename Lhs> 
inline const typename infernal::result_type<Lhs>::type 
    operator+(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(rhs) += lhs;    } 

// - Subtraction
template <typename Rhs> 
inline const typename infernal::result_type<Rhs>::type 
    operator-(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) -= rhs;    }

template <typename Lhs> 
inline const typename infernal::result_type<Lhs>::type 
    operator-(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(lhs) -= rhs;    }

// * Multiplication
template <typename Rhs> 
inline const typename infernal::result_type<Rhs>::type 
    operator*(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) *= rhs;    }

template <typename Lhs> 
inline const typename infernal::result_type<Lhs>::type 
    operator*(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(rhs) *= lhs;    } 

// / Division
template <typename Rhs> 
inline const typename infernal::result_type<Rhs>::type 
    operator/(const mpcomplex& lhs, const Rhs& rhs){ return mpcomplex(lhs) /= rhs;    }

template <typename Lhs> 
inline const typename infernal::result_type<Lhs>::type 
    operator/(const Lhs& lhs, const mpcomplex& rhs){ return mpcomplex(lhs) /= rhs;    }

//////////////////////////////////////////////////////////////////////////
// pow
const mpcomplex pow(const mpcomplex& a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const mpcomplex& a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const mpcomplex& a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const mpcomplex& a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const unsigned int a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const int a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long double a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const double a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const unsigned long int a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned long int a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned long int a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned long int a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned long int a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const unsigned int a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned int a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned int a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned int a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned int a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const unsigned int a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const long int a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long int a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const int a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const int a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const int a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const int a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd()); 
const mpcomplex pow(const int a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const int a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd()); 

const mpcomplex pow(const long double a, const long double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());    
const mpcomplex pow(const long double a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long double a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long double a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const long double a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

const mpcomplex pow(const double a, const double b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());    
const mpcomplex pow(const double a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const double a, const unsigned int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const double a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
const mpcomplex pow(const double a, const int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

inline const mpcomplex mul_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
inline const mpcomplex mul_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
inline const mpcomplex div_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());
inline const mpcomplex div_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd());

//////////////////////////////////////////////////////////////////////////
// 'Dirty' equality check 1: |a-b| < min{|a|,|b|} * eps
inline bool isEqualFuzzy(const mpcomplex& a, const mpcomplex& b, const mpreal& eps);

// 'Dirty' equality check 2: |a-b| < min{|a|,|b|} * eps( min{|a|,|b|} )
inline bool isEqualFuzzy(const mpcomplex& a, const mpcomplex& b);

// 'Dirty' realness check: |im(x)| < |x| * eps(|x|)
inline bool isRealFuzzy(const mpcomplex x);

// 'Bitwise' equality check
//  maxUlps - a and b can be apart by maxUlps binary numbers. 
inline bool isEqualUlps(const mpcomplex& a, const mpcomplex& b, int maxUlps);

//////////////////////////////////////////////////////////////////////////
// Convert precision in 'bits' to decimal digits and vice versa.
//    bits   = ceil(digits*log[2](10))
//    digits = floor(bits*log[10](2))

inline mp_prec_t digits2bits(int d);
inline int       bits2digits(mp_prec_t b);

//////////////////////////////////////////////////////////////////////////
// min, max
const mpcomplex (max)(const mpcomplex& x, const mpcomplex& y);
const mpcomplex (min)(const mpcomplex& x, const mpcomplex& y);

//////////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// Operators - Assignment
inline mpcomplex& mpcomplex::operator=(const mpcomplex& v)
{
    if (this != &v)
    {
		mp_prec_t tp = mpc_get_prec(  mpc_srcptr());
		mp_prec_t vp = mpc_get_prec(v.mpc_srcptr());

		if(tp != vp){
			clear(mpc_ptr());
			mpc_init2(mpc_ptr(), vp);
		}

        mpc_set(mpc_ptr(), v.mpc_srcptr(), mpcomplex::get_default_rnd());

        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    return *this;
}

template <typename real_t> 
inline mpcomplex& mpcomplex::operator= (const std::complex<real_t>& z)
{
	set_real(z.real());
	set_imag(z.imag());

	MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpreal& v)
{
	mp_prec_t tp = mpc_get_prec (  mpc_srcptr() );
	mp_prec_t vp = mpfr_get_prec(v.mpfr_srcptr());

	if(tp != vp){
		clear(mpc_ptr());
		mpc_init2(mpc_ptr(), vp);
	}

	mpc_set_fr(mpc_ptr(), v.mpfr_srcptr(), mpcomplex::get_default_rnd());

	MPREAL_MSVC_DEBUGVIEW_CODE;
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpf_t v)
{
    mpc_set_f(mpc_ptr(), v, mpcomplex::get_default_rnd());
    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const long double v)        
{    
    mpc_set_ld(mpc_ptr(), v, mpcomplex::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const double v)                
{   
#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
	if(fits_in_bits(v, MPREAL_DOUBLE_BITS_OVERFLOW))
	{
		mpc_set_d(mpc_ptr(),v,mpcomplex::get_default_rnd());
	}else
		throw conversion_overflow();
#else
	mpc_set_d(mpc_ptr(),v,mpcomplex::get_default_rnd());
#endif

	MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const unsigned long int v)    
{    
    mpc_set_ui(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const unsigned int v)        
{    
    mpc_set_ui(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const unsigned long long int v)    
{    
    mpc_set_uj(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const long long int v)    
{    
    mpc_set_sj(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const long int v)            
{    
    mpc_set_si(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const int v)
{    
    mpc_set_si(mpc_ptr(), v, mpcomplex::get_default_rnd());    

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const char* s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpcomplex(const char* s,        mp_prec_t prec, int base, mpc_rnd_t mode)
    //        mpcomplex(const std::string& s,mp_prec_t prec, int base, mpc_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpc_t t;

    mpc_init2(t, mpc_get_prec(mpc_srcptr()));

    if(0 == mpc_set_str(t, s, 10, mpcomplex::get_default_rnd()))
    {
        mpc_set(mpc_ptr(), t, mpcomplex::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const std::string& s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpcomplex(const char* s,        mp_prec_t prec, int base, mpc_rnd_t mode)
    //        mpcomplex(const std::string& s,mp_prec_t prec, int base, mpc_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpc_t t;

    mpc_init2(t, mpc_get_prec(mpc_srcptr()));

    if(0 == mpc_set_str(t, s.c_str(), 10, mpcomplex::get_default_rnd()))
    {
        mpc_set(mpc_ptr(), t, mpcomplex::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const mpreal& v)
{
	mpfr_set(mpc_re_ptr(),v.mpfr_srcptr(),mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const mpreal& v)
{
	mpfr_set(mpc_im_ptr(),v.mpfr_srcptr(),mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const mpf_t v)
{
	mpfr_set_f(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const mpf_t v)
{
	mpfr_set_f(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const long double v)
{
	mpfr_set_ld(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const long double v)
{
	mpfr_set_ld(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const double v)
{
	mpfr_set_d(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const double v)
{
	mpfr_set_d(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const unsigned long int v)
{
	mpfr_set_ui(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const unsigned long int v)
{
	mpfr_set_ui(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const unsigned long long int v)
{
	mpfr_set_uj(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const unsigned long long int v)
{
	mpfr_set_uj(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const long long int v)
{
	mpfr_set_sj(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const long long int v)
{
	mpfr_set_sj(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const unsigned int v)
{
	mpfr_set_ui(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const unsigned int v)
{
	mpfr_set_ui(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const long int v)
{
	mpfr_set_si(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const long int v)
{
	mpfr_set_si(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const int v)
{
	mpfr_set_si(mpc_re_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const int v)
{
	mpfr_set_si(mpc_im_ptr(),v,mpreal::get_default_rnd());

    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const char* s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s, 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpc_re_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
	mpfr_clear(t);
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const char* s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s, 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpc_im_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
	mpfr_clear(t);
    return *this;
}

inline mpcomplex& mpcomplex::set_real(const std::string& s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s.c_str(), 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpc_re_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
	mpfr_clear(t);
    return *this;
}

inline mpcomplex& mpcomplex::set_imag(const std::string& s)
{
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));
    if(0 == mpfr_set_str(t, s.c_str(), 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpc_im_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
	mpfr_clear(t);
    return *this;
}

//////////////////////////////////////////////////////////////////////////
// + Addition
inline mpcomplex& mpcomplex::operator+=(const mpcomplex& v)
{
    mpc_add(mpc_ptr(), mpc_srcptr(), v.mpc_srcptr(), mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpreal& v)
{
	mpc_add_fr(mpc_ptr(), mpc_srcptr(), v.mpfr_srcptr(), mpcomplex::get_default_rnd());
	MPREAL_MSVC_DEBUGVIEW_CODE;
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpf_t u)
{
    *this += mpcomplex(u);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+= (const long double u)
{
    *this += mpcomplex(u);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline mpcomplex& mpcomplex::operator+= (const double u)
{
    *this += mpcomplex(u);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const unsigned long int u)
{
    mpc_add_ui(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const unsigned int u)
{
    mpc_add_ui(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const long int u)
{
    mpc_add_si(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const int u)
{
    mpc_add_si(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const mpreal& v)
{
	mpfr_add(mpc_im_ptr(),mpc_im_srcptr(),v.mpfr_srcptr(),mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const mpf_t v)
{
	this->iplus(mpreal(v));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const long double u)
{
	mpfr_add_d(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const double u)
{
	mpfr_add_d(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const unsigned long int u)
{
	mpfr_add_ui(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const unsigned int u)
{
	mpfr_add_ui(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const long int u)
{
	mpfr_add_si(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iplus(const int u)
{
	mpfr_add_si(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator+=(const long long int u)         {    *this += mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::operator+=(const unsigned long long int u){    *this += mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::iplus(const long long int u)              {    this->iplus(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;   }
inline mpcomplex& mpcomplex::iplus(const unsigned long long int u)     {    this->iplus(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;   }
inline mpcomplex& mpcomplex::operator-=(const long long int  u)        {    *this -= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::operator-=(const unsigned long long int u){    *this -= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::iminus(const long long int u)             {    this->iminus(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;  }
inline mpcomplex& mpcomplex::iminus(const unsigned long long int u)    {    this->iminus(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;  }
inline mpcomplex& mpcomplex::operator*=(const long long int  u)        {    *this *= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::operator*=(const unsigned long long int u){    *this *= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::itimes(const long long int u)             {    this->itimes(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;  }
inline mpcomplex& mpcomplex::itimes(const unsigned long long int u)    {    this->itimes(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;  }
inline mpcomplex& mpcomplex::operator/=(const long long int  u)        {    *this /= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::operator/=(const unsigned long long int u){    *this /= mpcomplex(u); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::idiv(const long long int u)               {    this->idiv(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }
inline mpcomplex& mpcomplex::idiv(const unsigned long long int u)      {    this->idiv(mpreal(u)); MPREAL_MSVC_DEBUGVIEW_CODE; return *this;    }

inline const mpcomplex mpcomplex::operator+()const    {    return mpcomplex(*this); }

inline const mpcomplex operator+(const mpcomplex& a, const mpcomplex& b)
{
	mpcomplex c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
	mpc_add(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), mpcomplex::get_default_rnd());
	return c;
}

inline mpcomplex& mpcomplex::operator++() 
{
    return *this += 1;
}

inline const mpcomplex mpcomplex::operator++ (int)
{
    mpcomplex x(*this);
    *this += 1;
    return x;
}

inline mpcomplex& mpcomplex::operator--() 
{
    return *this -= 1;
}

inline const mpcomplex mpcomplex::operator-- (int)
{
    mpcomplex x(*this);
    *this -= 1;
    return x;
}

//////////////////////////////////////////////////////////////////////////
// - Subtraction
inline mpcomplex& mpcomplex::operator-=(const mpcomplex& v)
{
    mpc_sub(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator-=(const mpreal& v)
{
	mpc_sub_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),mpcomplex::get_default_rnd());
	MPREAL_MSVC_DEBUGVIEW_CODE;
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const long double v)
{
    *this -= mpcomplex(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline mpcomplex& mpcomplex::operator-=(const double v)
{
    *this -= mpcomplex(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator-=(const unsigned long int v)
{
    mpc_sub_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator-=(const unsigned int v)
{
    mpc_sub_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator-=(const long int v)
{
	mpreal V(v);
    mpc_sub_fr(mpc_ptr(),mpc_srcptr(),V.mpfr_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator-=(const int v)
{
	mpreal V(v);
    mpc_sub_fr(mpc_ptr(),mpc_srcptr(),V.mpfr_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const mpreal& v)
{
	mpfr_sub(mpc_im_ptr(),mpc_im_srcptr(),v.mpfr_srcptr(),mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const mpf_t v)
{
	this->iminus(mpreal(v));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const long double u)
{
	mpfr_sub_d(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const double u)
{
	mpfr_sub_d(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const unsigned long int u)
{
	mpfr_sub_ui(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const unsigned int u)
{
	mpfr_sub_ui(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const long int u)
{
	mpfr_sub_si(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::iminus(const int u)
{
	mpfr_sub_si(mpc_im_ptr(),mpc_im_srcptr(),u,mpreal::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const mpcomplex mpcomplex::operator-()const
{
    mpcomplex u(*this);
    mpc_neg(u.mpc_ptr(),u.mpc_srcptr(),mpcomplex::get_default_rnd());
    return u;
}

inline const mpcomplex operator-(const mpcomplex& a, const mpcomplex& b)
{
	mpcomplex c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
	mpc_sub(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), mpcomplex::get_default_rnd());
	return c;
}

inline const mpcomplex operator-(const double  b, const mpcomplex& a)
{
    mpcomplex x(b, mpc_get_prec(a.mpc_ptr()));
    x -= a;
    return x;
}

inline const mpcomplex operator-(const unsigned long int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_ptr()));
    mpc_ui_sub(x.mpc_ptr(), b, a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

inline const mpcomplex operator-(const unsigned int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_ptr()));
    mpc_ui_sub(x.mpc_ptr(), b, a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

inline const mpcomplex operator-(const long int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_ptr()));
	mpreal B(b);
    mpc_fr_sub(x.mpc_ptr(), B.mpfr_srcptr(), a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

inline const mpcomplex operator-(const int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_ptr()));
	mpreal B(b);
    mpc_fr_sub(x.mpc_ptr(), B.mpfr_srcptr(), a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

//////////////////////////////////////////////////////////////////////////
// * Multiplication
inline mpcomplex& mpcomplex::operator*= (const mpcomplex& v)
{
    mpc_mul(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator*=(const mpreal& v)
{
	mpc_mul_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),mpcomplex::get_default_rnd());
	MPREAL_MSVC_DEBUGVIEW_CODE;
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const long double v)
{
    *this *= mpcomplex(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline mpcomplex& mpcomplex::operator*=(const double v)
{
    *this *= mpcomplex(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator*=(const unsigned long int v)
{
    mpc_mul_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator*=(const unsigned int v)
{
    mpc_mul_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator*=(const long int v)
{
    mpc_mul_si(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator*=(const int v)
{
    mpc_mul_si(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const mpcomplex operator*(const mpcomplex& a, const mpcomplex& b)
{
	mpcomplex c(0, (std::max)(mpc_get_prec(a.mpc_ptr()), mpc_get_prec(b.mpc_ptr())));
	mpc_mul(c.mpc_ptr(), a.mpc_srcptr(), b.mpc_srcptr(), mpcomplex::get_default_rnd());
	return c;
}

inline mpcomplex& mpcomplex::itimes(const mpreal& v)
{
	mpc_mul_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const mpf_t v)
{
	this->itimes(mpreal(v));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const long double u)
{
	this->itimes(mpreal(u));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const double u)
{
	this->itimes(mpreal(u));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const unsigned long int u)
{
	mpc_mul_ui(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const unsigned int u)
{
	mpc_mul_ui(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const long int u)
{
	mpc_mul_si(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::itimes(const int u)
{
	mpc_mul_si(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

//////////////////////////////////////////////////////////////////////////
// / Division
inline mpcomplex& mpcomplex::operator/=(const mpcomplex& v)
{
    mpc_div(mpc_ptr(),mpc_srcptr(),v.mpc_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const mpreal& v)
{
    mpc_div_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const long double v)
{
    *this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;    
}

inline mpcomplex& mpcomplex::operator/=(const double v)
{
    *this /= mpreal(v);    
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const unsigned long int v)
{
    mpc_div_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const unsigned int v)
{
    mpc_div_ui(mpc_ptr(),mpc_srcptr(),v,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const long int v)
{
	*this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator/=(const int v)
{
	*this /= mpreal(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const mpreal& v)
{
	mpc_div_fr(mpc_ptr(),mpc_srcptr(),v.mpfr_srcptr(),get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),-1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const mpf_t v)
{
	this->idiv(mpreal(v));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const long double u)
{
	this->idiv(mpreal(u));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const double u)
{
	this->idiv(mpreal(u));
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const unsigned long int u)
{
	mpc_div_ui(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),-1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const unsigned int u)
{
	mpc_div_ui(mpc_ptr(),mpc_srcptr(),u,get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),-1,get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const long int u)
{
	mpc_div_ui(mpc_ptr(),mpc_srcptr(),std::abs(u),get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),(u<0)-(0<u),get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::idiv(const int u)
{
	mpc_div_ui(mpc_ptr(),mpc_srcptr(),std::abs(u),get_default_rnd());
	mpc_mul_i(mpc_ptr(), mpc_srcptr(),(u<0)-(0<u),get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const mpcomplex operator/(const mpcomplex& a, const mpcomplex& b)
{
	mpcomplex c(a, (std::max)(mpc_get_prec(a.mpc_srcptr()), mpc_get_prec(b.mpc_srcptr())));
	mpc_div(c.mpc_ptr(), c.mpc_srcptr(), b.mpc_srcptr(), mpcomplex::get_default_rnd());
	return c;
}

inline const mpcomplex operator/(const unsigned long int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_srcptr()));
    mpc_ui_div(x.mpc_ptr(), b, a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

inline const mpcomplex operator/(const unsigned int b, const mpcomplex& a)
{
    mpcomplex x(0, mpc_get_prec(a.mpc_srcptr()));
    mpc_ui_div(x.mpc_ptr(), b, a.mpc_srcptr(), mpcomplex::get_default_rnd());
    return x;
}

inline const mpcomplex operator/(const long int b, const mpcomplex& a)
{
    mpcomplex x(b, mpc_get_prec(a.mpc_srcptr()));
	x /= a;
    return x;
}

inline const mpcomplex operator/(const int b, const mpcomplex& a)
{
    mpcomplex x(b, mpc_get_prec(a.mpc_srcptr()));
	x /= a;
    return x;
}

inline const mpcomplex operator/(const double  b, const mpcomplex& a)
{
    mpcomplex x(b, mpc_get_prec(a.mpc_srcptr()));
    x /= a;
    return x;
}

//////////////////////////////////////////////////////////////////////////
// Shifts operators - Multiplication/Division by power of 2
inline mpcomplex& mpcomplex::operator<<=(const unsigned long int u)
{
    mpc_mul_2ui(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const unsigned int u)
{
    mpc_mul_2ui(mpc_ptr(),mpc_srcptr(),static_cast<unsigned long int>(u),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const long int u)
{
    mpc_mul_2si(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator<<=(const int u)
{
    mpc_mul_2si(mpc_ptr(),mpc_srcptr(),static_cast<long int>(u),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator>>=(const unsigned long int u)
{
    mpc_div_2ui(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator>>=(const unsigned int u)
{
    mpc_div_2ui(mpc_ptr(),mpc_srcptr(),static_cast<unsigned long int>(u),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator>>=(const long int u)
{
    mpc_div_2si(mpc_ptr(),mpc_srcptr(),u,mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::operator>>=(const int u)
{
    mpc_div_2si(mpc_ptr(),mpc_srcptr(),static_cast<long int>(u),mpcomplex::get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline const mpcomplex operator<<(const mpcomplex& v, const unsigned long int k)
{
    return mul_2ui(v,k);
}

inline const mpcomplex operator<<(const mpcomplex& v, const unsigned int k)
{
    return mul_2ui(v,static_cast<unsigned long int>(k));
}

inline const mpcomplex operator<<(const mpcomplex& v, const long int k)
{
    return mul_2si(v,k);
}

inline const mpcomplex operator<<(const mpcomplex& v, const int k)
{
    return mul_2si(v,static_cast<long int>(k));
}

inline const mpcomplex operator>>(const mpcomplex& v, const unsigned long int k)
{
    return div_2ui(v,k);
}

inline const mpcomplex operator>>(const mpcomplex& v, const long int k)
{
    return div_2si(v,k);
}

inline const mpcomplex operator>>(const mpcomplex& v, const unsigned int k)
{
    return div_2ui(v,static_cast<unsigned long int>(k));
}

inline const mpcomplex operator>>(const mpcomplex& v, const int k)
{
    return div_2si(v,static_cast<long int>(k));
}

// mul_2ui
inline const mpcomplex mul_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode)
{
    mpcomplex x(v);
    mpc_mul_2ui(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

// mul_2si
inline const mpcomplex mul_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode)
{
    mpcomplex x(v);
    mpc_mul_2si(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const mpcomplex div_2ui(const mpcomplex& v, unsigned long int k, mpc_rnd_t rnd_mode)
{
    mpcomplex x(v);
    mpc_div_2ui(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

inline const mpcomplex div_2si(const mpcomplex& v, long int k, mpc_rnd_t rnd_mode)
{
    mpcomplex x(v);
    mpc_div_2si(x.mpc_ptr(),v.mpc_srcptr(),k,rnd_mode);
    return x;
}

//////////////////////////////////////////////////////////////////////////
//Relational operators
inline bool operator == (const mpcomplex& a, const mpcomplex& b        ){  return (mpc_cmp(a.mpc_srcptr(),b.mpc_srcptr()) != 0 );                   }
inline bool operator == (const mpcomplex& a, const unsigned long int b ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const mpcomplex& a, const unsigned int b      ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const mpcomplex& a, const long int b          ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }
inline bool operator == (const mpcomplex& a, const int b               ){  return (mpc_cmp_si(a.mpc_srcptr(),b) == 0 );                }

inline bool operator != (const mpcomplex& a, const mpcomplex& b        ){  return !(a == b);  }
inline bool operator != (const mpcomplex& a, const unsigned long int b ){  return !(a == b);  }
inline bool operator != (const mpcomplex& a, const unsigned int b      ){  return !(a == b);  }
inline bool operator != (const mpcomplex& a, const long int b          ){  return !(a == b);  }
inline bool operator != (const mpcomplex& a, const int b               ){  return !(a == b);  }

//////////////////////////////////////////////////////////////////////////
// Type Converters
inline _Complex float mpcomplex::toCFloat (mpc_rnd_t mode) const
{
	_Complex float x = (_Complex float) mpc_get_dc(mpc_srcptr(), mode);
	return x;
}

inline std::complex<float> mpcomplex::toFloat (mpc_rnd_t mode) const
{
	_Complex float n = (_Complex float) mpc_get_dc(mpc_srcptr(), mode);
	std::complex<float> nxx = std::complex<float>(creal(n), cimag(n));
	return nxx;
}

inline _Complex double mpcomplex::toCDouble (mpc_rnd_t mode) const
{
	_Complex double x = mpc_get_dc(mpc_srcptr(), mode);
	return x;
}

inline std::complex<double> mpcomplex::toDouble (mpc_rnd_t mode) const
{
	_Complex double n = mpc_get_dc(mpc_srcptr(), mode);
	std::complex<double> nxx = std::complex<double>(creal(n), cimag(n));
	return nxx;
}

inline _Complex long double mpcomplex::toCLDouble (mpc_rnd_t mode) const
{
	_Complex long double x = mpc_get_ldc(mpc_srcptr(), mode);
	return x;
}

inline std::complex<long double> mpcomplex::toLDouble (mpc_rnd_t mode) const
{
	_Complex long double n = mpc_get_ldc(mpc_srcptr(), mode);
	std::complex<long double> nxx = std::complex<long double>(creal(n), cimag(n));
	return nxx;
}

inline ::mpc_ptr     mpcomplex::mpc_ptr()             { return mp; }
inline ::mpc_srcptr  mpcomplex::mpc_ptr()    const    { return mp; }
inline ::mpc_srcptr  mpcomplex::mpc_srcptr() const    { return mp; }

inline std::string mpcomplex::mpfrString(const ::mpfr_srcptr n, const std::string& format) const
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
inline std::string mpcomplex::toString(int n, int b, mpc_rnd_t mode) const 
{
	if(b == 10){
		std::string out("(");
    	std::ostringstream format;
	    int digits = (n >= 0) ? n : 1 + bits2digits(mpfr_get_prec(mpfr_srcptr()));
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

inline const mpcomplex mpc_const_infinity (int ReSign = 1, int ImSign = 0, mp_prec_t p = mpreal::get_default_prec())
{
	mpcomplex x(0, p);
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

//////////////////////////////////////////////////////////////////////////
// I/O
inline std::ostream& mpcomplex::output(std::ostream& os) const 
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

inline std::ostream& operator<<(std::ostream& os, const mpcomplex& v)
{
    return v.output(os);
}

inline std::istream& operator>>(std::istream &is, mpcomplex& v)
{
    // TODO: use cout::hexfloat and other flags to setup base
    std::string tmp;
    is >> tmp;
    mpc_set_str(v.mpc_ptr(), tmp.c_str(), 10, mpcomplex::get_default_rnd());
    return is;
}

//////////////////////////////////////////////////////////////////////////
// Set/Get number properties
inline int mpcomplex::getPrecision() const
{
    return int(mpc_get_prec(mpc_srcptr()));
}

inline mpcomplex& mpcomplex::setPrecision(int Precision, mpc_rnd_t RoundingMode)
{
	mpc_t v;
	mpc_init2(v, Precision);
	mpc_set(v, mpc_srcptr(), RoundingMode);
	mpc_swap(v, mpc_ptr());
	clear(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mp_prec_t mpcomplex::get_prec() const
{
    return mpc_get_prec(mpc_srcptr());
}

inline void mpcomplex::set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode)
{
	mpc_t v;
	mpc_init2(v, prec);
	mpc_set(v, mpc_srcptr(), rnd_mode);
	mpc_swap(v, mpc_ptr());
	clear(v);
    MPREAL_MSVC_DEBUGVIEW_CODE;
}

inline mpcomplex& mpcomplex::setInf(int ReSign, int ImSign)
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

inline mpcomplex& mpcomplex::setNan() 
{
    mpc_set_nan(mpc_ptr());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::setZero(int ReSign, int ImSign)
{
    mpc_set_si(mpc_ptr(), 0, get_default_rnd());
    setSigns(ReSign, ImSign, get_default_rnd());
    MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

inline mpcomplex& mpcomplex::setSigns(int ReSign, int ImSign, mpc_rnd_t rnd_mode)
{
	mpfr_setsign(mpc_re_ptr(), mpc_re_srcptr(), ReSign, MPC_RND_RE(rnd_mode));
	mpfr_setsign(mpc_im_ptr(), mpc_im_srcptr(), ImSign, MPC_RND_IM(rnd_mode));
    MPREAL_MSVC_DEBUGVIEW_CODE;
	return *this;
}

inline bool isEqualUlps(const mpcomplex& a, const mpcomplex& b, int maxUlps)
{
    return abs(a - b) <= machine_epsilon((max)(abs(a), abs(b))) * maxUlps;
}

inline bool isEqualFuzzy(const mpcomplex& a, const mpcomplex& b, const mpreal& eps)
{
    return abs(a - b) <= eps;
}

inline bool isEqualFuzzy(const mpcomplex& a, const mpcomplex& b)
{
    return isEqualFuzzy(a, b, machine_epsilon((max)(1, (min)(abs(a), abs(b)))));
}

inline bool isRealFuzzy(const mpcomplex x)
{
	return isEqualFuzzy(imag(x), mpreal(0), machine_epsilon(abs(x)));
}

//////////////////////////////////////////////////////////////////////////
// Mathematical Functions
//////////////////////////////////////////////////////////////////////////
#define MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(f)                    \
        mpcomplex y(0, mpc_get_prec(x.mpc_srcptr()));          \
        mpc_##f(y.mpc_ptr(), x.mpc_srcptr(), r);           \
        return y; 

inline const mpcomplex sqr  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd())
{   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqr );    }

inline const mpcomplex sqrt (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd())
{   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqrt);    }

// ~this probably shouldn't be here since it's not an mpc base function
inline const mpcomplex root(const mpcomplex& x, unsigned long int k, mpc_rnd_t r = mpcomplex::get_default_rnd())
{
    mpcomplex y(0, mpc_get_prec(x.mpc_srcptr())); 
    mpreal inv(k, mpc_get_prec(x.mpc_srcptr()), MPC_RND_RE(r));
    inv = 1/inv;
    mpc_pow_fr(y.mpc_ptr(), x.mpc_srcptr(), inv.mpfr_srcptr(), r);
    return y; 
}

inline int sin_cos(mpcomplex& s, mpcomplex& c, const mpcomplex& v, mpc_rnd_t rnd_mode_s = mpcomplex::get_default_rnd(), mpc_rnd_t rnd_mode_c = mpcomplex::get_default_rnd())
{
    return mpc_sin_cos(s.mpc_ptr(), c.mpc_ptr(), v.mpc_srcptr(), rnd_mode_s, rnd_mode_c);
}

inline const mpreal abs (const mpcomplex& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_abs(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal norm (const mpcomplex& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_norm(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal real (const mpcomplex& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_real(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal imag (const mpcomplex& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_imag(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpreal arg  (const mpcomplex& x, mpc_rnd_t r)
{
    mpreal y(0, x.get_prec(), MPC_RND_RE(r));
    mpc_arg(y.mpfr_ptr(), x.mpc_srcptr(), MPC_RND_RE(r));
    return y;
}

inline const mpcomplex conj (const mpcomplex& x, mpc_rnd_t r)
{
    mpcomplex y(x);
    mpc_conj(y.mpc_ptr(), x.mpc_srcptr(), r);
    return y;
}

inline const mpcomplex sqrt  (const mpreal v, mpc_rnd_t rnd_mode)               {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const unsigned long int v, mpc_rnd_t rnd_mode)    {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const unsigned int v, mpc_rnd_t rnd_mode)         {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const long int v, mpc_rnd_t rnd_mode)             {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const int v, mpc_rnd_t rnd_mode)                  {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const long double v, mpc_rnd_t rnd_mode)          {   return sqrt(mpcomplex(v),rnd_mode);    }
inline const mpcomplex sqrt  (const double v, mpc_rnd_t rnd_mode)               {   return sqrt(mpcomplex(v),rnd_mode);    }

inline const mpcomplex log   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(log  );    }
inline const mpcomplex log10 (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(log10);    }
inline const mpcomplex exp   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(exp  );    }
inline const mpcomplex cos   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cos  );    }
inline const mpcomplex sin   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sin  );    }
inline const mpcomplex tan   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tan  );    }
inline const mpcomplex acos  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acos );    }
inline const mpcomplex asin  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asin );    }
inline const mpcomplex atan  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atan );    }
inline const mpcomplex cosh  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cosh );    }
inline const mpcomplex sinh  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sinh );    }
inline const mpcomplex tanh  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tanh );    }
inline const mpcomplex acosh (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acosh);    }
inline const mpcomplex asinh (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asinh);    }
inline const mpcomplex atanh (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atanh);    }

inline const mpcomplex sec   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/cos(x, r);                           }
inline const mpcomplex csc   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/sin(x, r);                           }
inline const mpcomplex cot   (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/tan(x, r);                           }
inline const mpcomplex sech  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/cosh(x, r);                          }
inline const mpcomplex csch  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/sinh(x, r);                          }
inline const mpcomplex coth  (const mpcomplex& x, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return 1/tanh(x, r);                          }
inline const mpcomplex acot  (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return atan (1/v, r);                      }
inline const mpcomplex asec  (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return acos (1/v, r);                      }
inline const mpcomplex acsc  (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return asin (1/v, r);                      }
inline const mpcomplex acoth (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return atanh(1/v, r);                      }
inline const mpcomplex asech (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return acosh(1/v, r);                      }
inline const mpcomplex acsch (const mpcomplex& v, mpc_rnd_t r = mpcomplex::get_default_rnd()) {   return asinh(1/v, r);                      }

inline const mpcomplex fma (const mpcomplex& v1, const mpcomplex& v2, const mpcomplex& v3, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd())
{
    mpcomplex a;
    mp_prec_t p1, p2, p3;

    p1 = v1.get_prec(); 
    p2 = v2.get_prec(); 
    p3 = v3.get_prec(); 

    a.set_prec(p3>p2?(p3>p1?p3:p1):(p2>p1?p2:p1));

    mpc_fma(a.mpc_ptr(),v1.mpc_srcptr(),v2.mpc_srcptr(),v3.mpc_srcptr(),rnd_mode);
    return a;
}

//////////////////////////////////////////////////////////////////////////
// Miscellaneous Functions
inline void         swap (mpcomplex& a, mpcomplex& b)            {    mpc_swap(a.mp,b.mp);   }

inline const mpcomplex mpc_urandom (gmp_randstate_t& state)
{
    mpcomplex x;
    mpc_urandom(x.mpc_ptr(), state);
    return x;
}

//////////////////////////////////////////////////////////////////////////
// Set/Get global properties
inline void mpcomplex::set_default_prec(mp_prec_t prec)
{ 
    mpfr_set_default_prec(prec); 
}

inline bool mpcomplex::fits_in_bits(double x, int n)
{   
    int i;
    double t;
    return IsInf(x) || (std::modf ( std::ldexp ( std::frexp ( x, &i ), n ), &t ) == 0.0);
}

inline const mpcomplex pow(const mpcomplex& a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd())
{
    mpcomplex x(a);
    mpc_pow(x.mp,x.mp,b.mp,rnd_mode);
    return x;
}

inline const mpcomplex pow(const mpcomplex& a, const unsigned long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd())
{
    mpcomplex x(a);
    mpc_pow_ui(x.mp,x.mp,b,rnd_mode);
    return x;
}

inline const mpcomplex pow(const mpcomplex& a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(a,static_cast<unsigned long int>(b),rnd_mode);
}

inline const mpcomplex pow(const mpcomplex& a, const long int b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd())
{
    mpcomplex x(a);
    mpc_pow_si(x.mp,x.mp,b,rnd_mode);
    return x;
}

inline const mpcomplex pow(const mpcomplex& a, const int b, mpc_rnd_t rnd_mode)
{
    return pow(a,static_cast<long int>(b),rnd_mode);
}

inline const mpcomplex pow(const mpcomplex& a, const long double b, mpc_rnd_t rnd_mode)
{
    return pow(a,mpcomplex(b),rnd_mode);
}

inline const mpcomplex pow(const mpcomplex& a, const double b, mpc_rnd_t rnd_mode)
{
    return pow(a,mpcomplex(b),rnd_mode);
}

inline const mpcomplex pow(const unsigned long int a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::get_default_rnd())
{
    mpcomplex x(a);
    mpc_pow(x.mpc_ptr(),x.mpc_srcptr(),b.mpc_srcptr(),rnd_mode);
    return x;
}

inline const mpcomplex pow(const unsigned int a, const mpcomplex& b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),b,rnd_mode);
}

inline const mpcomplex pow(const long int a, const mpcomplex& b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),b,rnd_mode);
}

inline const mpcomplex pow(const int a, const mpcomplex& b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),b,rnd_mode);
}

inline const mpcomplex pow(const long double a, const mpcomplex& b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode);
}

inline const mpcomplex pow(const double a, const mpcomplex& b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode);
}

// pow unsigned long int
inline const mpcomplex pow(const unsigned long int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    mpreal x(a);
    mpfr_ui_pow_ui(x.mpfr_ptr(),a,b,MPC_RND_RE(rnd_mode));
    return mpcomplex(x);
}

inline const mpcomplex pow(const unsigned long int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const mpcomplex pow(const unsigned long int a, const long int b, mpc_rnd_t rnd_mode)
{
    if(b>0)    return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else       return pow(a,mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned long int a, const int b, mpc_rnd_t rnd_mode)
{
    if(b>0)    return pow(a,static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else       return pow(a,mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned long int a, const long double b, mpc_rnd_t rnd_mode)
{
    return pow(a,mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned long int a, const double b, mpc_rnd_t rnd_mode)
{
    return pow(a,mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

// pow unsigned int
inline const mpcomplex pow(const unsigned int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
}

inline const mpcomplex pow(const unsigned int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
}

inline const mpcomplex pow(const unsigned int a, const long int b, mpc_rnd_t rnd_mode)
{
    if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned int a, const int b, mpc_rnd_t rnd_mode)
{
    if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
    else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned int a, const long double b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

inline const mpcomplex pow(const unsigned int a, const double b, mpc_rnd_t rnd_mode)
{
    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
}

// pow long int
inline const mpcomplex pow(const long int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    if (a>0) return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
    else     return pow(mpcomplex(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const long int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    if (a>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
    else     return pow(mpcomplex(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const long int a, const long int b, mpc_rnd_t rnd_mode)
{
    if (a>0)
    {
        if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    }else{
        return pow(mpcomplex(a),b,rnd_mode); // mpfr_pow_si
    }
}

inline const mpcomplex pow(const long int a, const int b, mpc_rnd_t rnd_mode)
{
    if (a>0)
    {
        if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    }else{
        return pow(mpcomplex(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
    }
}

inline const mpcomplex pow(const long int a, const long double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    else        return pow(mpcomplex(a),mpcomplex(b),rnd_mode); //mpfr_pow
}

inline const mpcomplex pow(const long int a, const double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    else        return pow(mpcomplex(a),mpcomplex(b),rnd_mode); //mpfr_pow
}

// pow int
inline const mpcomplex pow(const int a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    if (a>0) return pow(static_cast<unsigned long int>(a),b,rnd_mode); //mpfr_ui_pow_ui
    else     return pow(mpcomplex(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const int a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    if (a>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode);  //mpfr_ui_pow_ui
    else     return pow(mpcomplex(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const int a, const long int b, mpc_rnd_t rnd_mode)
{
    if (a>0)
    {
        if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    }else{
        return pow(mpcomplex(a),b,rnd_mode); // mpfr_pow_si
    }
}

inline const mpcomplex pow(const int a, const int b, mpc_rnd_t rnd_mode)
{
    if (a>0)
    {
        if(b>0) return pow(static_cast<unsigned long int>(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_ui_pow_ui
        else    return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    }else{
        return pow(mpcomplex(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
    }
}

inline const mpcomplex pow(const int a, const long double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    else        return pow(mpcomplex(a),mpcomplex(b),rnd_mode); //mpfr_pow
}

inline const mpcomplex pow(const int a, const double b, mpc_rnd_t rnd_mode)
{
    if (a>=0)   return pow(static_cast<unsigned long int>(a),mpcomplex(b),rnd_mode); //mpfr_ui_pow
    else        return pow(mpcomplex(a),mpcomplex(b),rnd_mode); //mpfr_pow
}

// pow long double 
inline const mpcomplex pow(const long double a, const long double b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),mpcomplex(b),rnd_mode);
}

inline const mpcomplex pow(const long double a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const long double a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),static_cast<unsigned long int>(b),rnd_mode); //mpfr_pow_ui
}

inline const mpcomplex pow(const long double a, const long int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode); // mpfr_pow_si
}

inline const mpcomplex pow(const long double a, const int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}

inline const mpcomplex pow(const double a, const double b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),mpcomplex(b),rnd_mode);
}

inline const mpcomplex pow(const double a, const unsigned long int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode); // mpfr_pow_ui
}

inline const mpcomplex pow(const double a, const unsigned int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),static_cast<unsigned long int>(b),rnd_mode); // mpfr_pow_ui
}

inline const mpcomplex pow(const double a, const long int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),b,rnd_mode); // mpfr_pow_si
}

inline const mpcomplex pow(const double a, const int b, mpc_rnd_t rnd_mode)
{
    return pow(mpcomplex(a),static_cast<long int>(b),rnd_mode); // mpfr_pow_si
}
} // End of mpfr namespace

// Explicit specialization of std::swap for mpcomplex numbers
// Thus standard algorithms will use efficient version of swap (due to Koenig lookup)
// Non-throwing swap C++ idiom: http://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Non-throwing_swap
namespace std
{
	// we are allowed to extend namespace std with specializations only
    template <>
    inline void swap(mpfr::mpcomplex& x, mpfr::mpcomplex& y) 
    { 
        return mpfr::swap(x, y); 
    }

    template<>
    class numeric_limits<mpfr::mpcomplex>
    {
    public:
        static const bool is_specialized    = true;
        static const bool is_signed         = true;
        static const bool is_integer        = false;
        static const bool is_exact          = false;
        static const int  radix             = 2;    

        static const bool has_infinity      = true;
        static const bool has_quiet_NaN     = true;
        static const bool has_signaling_NaN = true;

        static const bool is_iec559         = true;        // = IEEE 754
        static const bool is_bounded        = true;
        static const bool is_modulo         = false;
        static const bool traps             = true;
        static const bool tinyness_before   = true;

        static const float_denorm_style has_denorm  = denorm_absent;

        // Returns smallest eps such that x + eps != x (relative machine epsilon)
		// NOTE: THIS SHOULD USE MAX OF RESULTS FOR REAL AND IMAG PARTS
        inline static mpfr::mpcomplex epsilon(const mpfr::mpcomplex& x) {  return mpfr::machine_epsilon(abs(x));  }

        inline static mpfr::mpcomplex round_error(mp_prec_t precision = mpfr::mpcomplex::get_default_prec())
        {
            mp_rnd_t r = mpfr::mpreal::get_default_rnd();

            if(r == GMP_RNDN)  return mpfr::mpcomplex(0.5, precision); 
            else               return mpfr::mpcomplex(1.0, precision);    
        }

        inline static const mpfr::mpcomplex infinity()         { return mpfr::mpc_const_infinity();     }
        inline static const mpfr::mpcomplex quiet_NaN()        { return mpfr::mpcomplex().setNan();    }
        inline static const mpfr::mpcomplex signaling_NaN()    { return mpfr::mpcomplex().setNan();    }

        // Please note, exponent range is not fixed in MPFR
        static const int min_exponent = MPFR_EMIN_DEFAULT;
        static const int max_exponent = MPFR_EMAX_DEFAULT;
        MPREAL_PERMISSIVE_EXPR static const int min_exponent10 = (int) (MPFR_EMIN_DEFAULT * 0.3010299956639811); 
        MPREAL_PERMISSIVE_EXPR static const int max_exponent10 = (int) (MPFR_EMAX_DEFAULT * 0.3010299956639811); 

#ifdef MPREAL_HAVE_DYNAMIC_STD_NUMERIC_LIMITS

        // Following members should be constant according to standard, but they can be variable in MPFR
        // So we define them as functions here. 
        //
        // This is preferable way for std::numeric_limits<mpfr::mpcomplex> specialization.
        // But it is incompatible with standard std::numeric_limits and might not work with other libraries, e.g. boost. 
        // See below for compatible implementation. 
        inline static float_round_style round_style()
        {
            mp_rnd_t r = mpfr::mpreal::get_default_rnd();

            switch (r)
            {
            case GMP_RNDN: return round_to_nearest;
            case GMP_RNDZ: return round_toward_zero; 
            case GMP_RNDU: return round_toward_infinity; 
            case GMP_RNDD: return round_toward_neg_infinity; 
            default: return round_indeterminate;
            }
        }

        inline static int digits()                        {    return int(mpfr::mpcomplex::get_default_prec());    }
        inline static int digits(const mpfr::mpcomplex& x)   {    return x.getPrecision();                         }

        inline static int digits10(mp_prec_t precision = mpfr::mpcomplex::get_default_prec())
        {
            return mpfr::bits2digits(precision);
        }

        inline static int digits10(const mpfr::mpcomplex& x)
        {
            return mpfr::bits2digits(x.getPrecision());
        }

        inline static int max_digits10(mp_prec_t precision = mpfr::mpcomplex::get_default_prec())
        {
            return digits10(precision);
        }
#else
        // Digits and round_style are NOT constants when it comes to mpcomplex.
        // If possible, please use functions digits() and round_style() defined above.
        //
        // These (default) values are preserved for compatibility with existing libraries, e.g. boost.
        // Change them accordingly to your application. 
        //
        // For example, if you use 256 bits of precision uniformly in your program, then:
        // digits       = 256
        // digits10     = 77 
        // max_digits10 = 78
        // 
        // Approximate formula for decimal digits is: digits10 = floor(log10(2) * digits). See bits2digits() for more details.

        static const std::float_round_style round_style = round_to_nearest;
        static const int digits       = 768;
        static const int digits10     = 231;
        static const int max_digits10 = 232;
#endif
    };

}

#endif /* __MPCOMPLEX_H__ */
