/*
** cmath.c - Math module with complex numbers
**
** See Copyright Notice in mruby.h
*/

/*
** This `mruby-cmath-alt` gem uses C99 _Complex features and GCC extensions,
** but not complex.h.
** You need a version of GCC that supports C99+.
*/

#include <mruby.h>

#ifdef MRB_NO_FLOAT
# error CMath conflicts with 'MRB_NO_FLOAT' configuration
#endif

mrb_value mrb_complex_new(mrb_state *mrb, mrb_float real, mrb_float imag);
void mrb_complex_get(mrb_state *mrb, mrb_value cpx, mrb_float*, mrb_float*);

static mrb_bool
cmath_get_complex(mrb_state *mrb, mrb_value c, mrb_float *r, mrb_float *i)
{
  if (mrb_integer_p(c)) {
    *r = (mrb_float)mrb_integer(c);
    *i = 0;
    return FALSE;
  }
  else if (mrb_float_p(c)) {
    *r = mrb_float(c);
    *i = 0;
    return FALSE;
  }
  else if (mrb_type(c) == MRB_TT_COMPLEX) {
    mrb_complex_get(mrb, c, r, i);
    return TRUE;
  }
  else {
    mrb_raise(mrb, E_TYPE_ERROR, "Numeric required");
    return FALSE;
  }
}

#ifdef MRB_USE_FLOAT32
#define F(x) x##f
#else
#define F(x) x
#endif

#if defined(_WIN32) && !defined(__MINGW32__)

#ifdef MRB_USE_FLOAT32
typedef _Fcomplex mrb_complex;
#else
typedef _Dcomplex mrb_complex;
#endif

static mrb_complex
CXDIVf(mrb_complex x, mrb_float y)
{
  return cmath_build_complex(cmath_creal(x)/y, cmath_cimag(x)/y);
}

static mrb_complex
CXDIVc(mrb_complex a, mrb_complex b)
{
  mrb_float ratio, den;
  mrb_float abr, abi, cr, ci;

  if ((abr = cmath_creal(b)) < 0)
    abr = - abr;
  if ((abi = cmath_cimag(b)) < 0)
    abi = - abi;
  if (abr <= abi) {
    ratio = cmath_creal(b) / cmath_cimag(b);
    den = cmath_cimag(a) * (1 + ratio*ratio);
    cr = (cmath_creal(a)*ratio + cmath_cimag(a)) / den;
    ci = (cmath_cimag(a)*ratio - cmath_creal(a)) / den;
  }
  else {
    ratio = cmath_cimag(b) / cmath_creal(b);
    den = cmath_creal(a) * (1 + ratio*ratio);
    cr = (cmath_creal(a) + cmath_cimag(a)*ratio) / den;
    ci = (cmath_cimag(a) - cmath_creal(a)*ratio) / den;
  }
  return cmath_build_complex(cr, ci);
}

#else

#if defined(__cplusplus) && (defined(__APPLE__) || (defined(__clang__) && (defined(__FreeBSD__) || defined(__OpenBSD__))))

#ifdef MRB_USE_FLOAT32
typedef std::complex<float> mrb_complex;
#else
typedef std::complex<double> mrb_complex;
#endif  /* MRB_USE_FLOAT32 */


#else  /* cpp */

#ifdef MRB_USE_FLOAT32
typedef float _Complex mrb_complex;
#else
typedef double _Complex mrb_complex;
#endif  /*  MRB_USE_FLOAT32 */

static mrb_complex
cmath_build_complex(mrb_float x, mrb_float y)
{
#ifdef __GNUC__
  return __builtin_complex(x, y);
#else
  union { mrb_float r[2]; mrb_complex c; } u;

  u.r[0] = x;
  u.r[1] = y;
  return u.c;
#endif
}

static mrb_float
cmath_creal(mrb_complex c)
{
#ifdef __GNUC__
  return __real__(c);
#else
  union { mrb_float r[2]; mrb_complex c; } u;

  u.c = c;
  return u.r[0];
#endif
}

static mrb_float
cmath_cimag(mrb_complex c)
{
#ifdef __GNUC__
  return __imag__(c);
#else
  union { mrb_float r[2]; mrb_complex c; } u;

  u.c = c;
  return u.r[1];
#endif
}
#endif

#define CXDIVf(x,y) (x)/(y)
#define CXDIVc(x,y) (x)/(y)

#endif

#define DEF_CMATH_METHOD(name) \
static mrb_value \
cmath_ ## name(mrb_state *mrb, mrb_value self)\
{\
  mrb_value z = mrb_get_arg1(mrb);\
  mrb_float real, imag;\
  if (cmath_get_complex(mrb, z, &real, &imag)) {\
    mrb_complex c = cmath_build_complex(real,imag);\
    c = cmath_c ## name(c);\
    return mrb_complex_new(mrb, cmath_creal(c), cmath_cimag(c));\
  }\
  return mrb_float_value(mrb, F(name)(real));\
}

static mrb_complex
cmath_cexp(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);

  if (isnan(x)) {
    if (y == 0.0F) {
      return cmath_build_complex(NAN, y);
    } else {
      return cmath_build_complex(NAN, NAN);
    }
  }
  if (x == +INFINITY) {
    if (isnan(y) || isinf(y)) {
      return cmath_build_complex(+INFINITY, NAN);
    } else if (y == 0.0F) {
      return c;
    }
  } else if (x == -INFINITY) {
    if (isnan(y) || isinf(y)) {
      return cmath_build_complex(+0.0F, F(copysign)(0.0F, y));
    }
  }

  mrb_float r = F(exp)(x);
  return cmath_build_complex(r*F(cos)(y), r*F(sin)(y));
}

static mrb_complex
cmath_clog(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float r = F(hypot)(x, y);
  mrb_float t = F(atan2)(y, x);
  return cmath_build_complex(F(log)(r), t);
}

static mrb_complex
cmath_csqrt(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);

  if (y == 0.0F) {
    if (isnan(x)) {
      return cmath_build_complex(x, x);
    } else if (signbit(x)) {
      return cmath_build_complex(0.0F, F(copysign)(F(sqrt)(-x), y));
    } else {
      return cmath_build_complex(F(sqrt)(+x), y);
    }
  } else {
    if (isinf(x) && isinf(y)) {
      return cmath_build_complex(INFINITY, y);
    } else if (isinf(x) && isnan(y)) {
      if (signbit(x)) {
        return cmath_build_complex(y, INFINITY);
      } else {
        return c;
      }
    } else if (isinf(x)) {
      if (signbit(x)) {
        return cmath_build_complex(0.0F, F(copysign)(INFINITY, y));
      } else {
        return cmath_build_complex(INFINITY, F(copysign)(0.0F, y));
      }
    } else if (isinf(y)) {
      return cmath_build_complex(INFINITY, y);
    } else {
#ifdef MRB_USE_FLOAT32
      static const float cutoff = 1e38;
#else
      static const double cutoff = 1e308;
#endif
      _Bool scale = (F(fabs)(x) > cutoff || (F(fabs)(y) > cutoff));
      if (scale) {
        /* Prevent hypot from overflowing */
        x /= 4.0F;
        y /= 4.0F;
      }
      mrb_float r = F(hypot)(x, y);
      mrb_float t = F(atan2)(y, x);
      r = F(sqrt)(r);
      t /= 2.0F;
      if (scale) {
        r *= 2.0F;
      }
      return cmath_build_complex(r*F(cos)(t), r*F(sin)(t));
    }
  }
}

static mrb_complex
cmath_csin(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float cx = F(cos)(x);
  mrb_float sx = F(sin)(x);
  mrb_float cy = F(cosh)(y);
  mrb_float sy = F(sinh)(y);
  return cmath_build_complex(sx*cy, cx*sy);
}

static mrb_complex
cmath_ccos(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float cx = F(cos)(x);
  mrb_float sx = F(sin)(x);
  mrb_float cy = F(cosh)(y);
  mrb_float sy = F(sinh)(y);
  return cmath_build_complex(cx*cy, -sx*sy);
}

static mrb_complex
cmath_ctan(mrb_complex c)
{
#ifdef MRB_USE_FLOAT32
  static const float cutoff1 = 53.0F;
  static const float cutoff2 = 0x1.0A2B24P+3F;
#else
  static const double cutoff1 = 373.0;
  static const double cutoff2 = 0x1.3001004048044P+4;
#endif
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float cx = F(cos)(x);
  mrb_float sx = F(sin)(x);
  mrb_complex w;

  if (F(fabs)(y) > cutoff1) {
    /* Cutoff above which real(w) == 0.0 */
    w = cmath_build_complex(F(copysign)(0.0F, sx*cx), F(copysign)(1.0F, y));
  } else if (F(fabs)(y) > cutoff2) {
    /* Cutoff above which |sy| == cy */
    mrb_float cy = F(cosh)(y);
    /* Not (sx*cx)/(cy*cy); cy*cy might overflow */
    w = cmath_build_complex(sx*cx/cy/cy, F(copysign)(1.0F, y));
  } else {
    mrb_float cy = F(cosh)(y);
    mrb_float sy = F(sinh)(y);
    mrb_float d = cx*cx*cy*cy + sx*sx*sy*sy;
    w = cmath_build_complex(sx*cx/d, sy*cy/d);
  }
  return w;
}

static mrb_complex
cmath_csinh(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float cx = F(cosh)(x);
  mrb_float sx = F(sinh)(x);
  mrb_float cy = F(cos)(y);
  mrb_float sy = F(sin)(y);
  return cmath_build_complex(sx*cy, cx*sy);
}

static mrb_complex
cmath_ccosh(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  if (isnan(x)) {
    if (isnan(y) || isinf(y)) {
      return cmath_build_complex(NAN, NAN);
    } else {
      return cmath_build_complex(NAN, y == 0.0F ? y : NAN);
    }
  } else if (isinf(x)) {
    if (isnan(y) || isinf(y)) {
      return cmath_build_complex(INFINITY, NAN);
    } else if (y == 0.0F) {
      return cmath_build_complex(INFINITY, signbit(x) ? -y : +y);
    } else {
      mrb_float cy = F(cos)(y);
      mrb_float sy = F(sin)(y);
      return cmath_build_complex(INFINITY*cy, x*sy);
    }
  } else {
    if (isnan(y) || isinf(y)) {
      return cmath_build_complex(NAN, x == 0.0F ? 0.0F : NAN);
    } else {
      mrb_float cx = F(cosh)(x);
      mrb_float sx = F(sinh)(x);
      mrb_float cy = F(cos)(y);
      mrb_float sy = F(sin)(y);
      return cmath_build_complex(cx*cy, sx*sy);
    }
  }
}

static mrb_complex
cmath_ctanh(mrb_complex c)
{
#ifdef MRB_USE_FLOAT32
  static const float cutoff1 = 53.0F;
  static const float cutoff2 = 0x1.0A2B24P+3F;
#else
  static const double cutoff1 = 373.0;
  static const double cutoff2 = 0x1.3001004048044P+4;
#endif
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);
  mrb_float cy = F(cos)(y);
  mrb_float sy = F(sin)(y);
  mrb_complex w;

  if (F(fabs)(x) > cutoff1) {
    /* Cutoff above which imag(w) == 0.0 */
    w = cmath_build_complex(F(copysign)(1.0F, x), 0.0F);
  } else if (F(fabs)(x) > cutoff2) {
    /* Cutoff above which |sx| == cx */
    mrb_float cx = F(cosh)(x);
    /* Not (sy*cy)/(cx*cx); cx*cx might overflow */
    w = cmath_build_complex(F(copysign)(1.0F, x), sy*cy/cx/cx);
  } else {
    mrb_float cx = F(cosh)(x);
    mrb_float sx = F(sinh)(x);
    mrb_float d = cx*cx*cy*cy + sx*sx*sy*sy;
    w = cmath_build_complex(sx*cx/d, sy*cy/d);
  }
  return w;
}

static mrb_complex
cmath_casinh(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);

  if (F(fabs)(x) > 1e8F || F(fabs)(y) > 1e8F) {
    /* Above this cutoff, c*c+1 == c*c; below it, c*c never overflows */
    if (signbit(x)) {
      return -(cmath_clog(-c) + (mrb_float)0.69314718055994530942);
    } else {
      return +(cmath_clog(+c) + (mrb_float)0.69314718055994530942);
    }
  } else {
    return cmath_clog(c + cmath_csqrt(c*c + 1.0F));
  }
}

static mrb_complex
cmath_cacosh(mrb_complex c)
{
  mrb_float x = cmath_creal(c);
  mrb_float y = cmath_cimag(c);

  if (F(fabs)(x) > 1e8F || F(fabs)(y) > 1e8F) {
    /* Above this cutoff, c*c-1 == c*c; below it, c*c never overflows */
    return cmath_clog(c) + (mrb_float)0.69314718055994530942;
  } else {
    return cmath_clog(c + cmath_csqrt(c + 1.0F)*cmath_csqrt(c - 1.0F));
  }
}

static mrb_complex
cmath_catanh(mrb_complex c)
{
  return 0.5F*cmath_clog((1.0F + c)/(1.0F - c));
}

static mrb_complex
cmath_casin(mrb_complex c)
{
  /* -i*asinh(i*c) */
  mrb_float x1 = cmath_creal(c);
  mrb_float y1 = cmath_cimag(c);
  mrb_complex c2 = cmath_build_complex(-y1, +x1);
  mrb_complex d2 = cmath_casinh(c2);
  mrb_float x2 = cmath_creal(d2);
  mrb_float y2 = cmath_cimag(d2);
  return cmath_build_complex(+y2, -x2);
}

static mrb_complex
cmath_cacos(mrb_complex c)
{
  /* -i*acosh(c) */
  mrb_complex d2 = cmath_cacosh(c);
  mrb_float x2 = cmath_creal(d2);
  mrb_float y2 = cmath_cimag(d2);
  mrb_complex d = cmath_build_complex(+y2, -x2);
  if (signbit(cmath_creal(d))) {
    d = -d;
  }
  return d;
}

static mrb_complex
cmath_catan(mrb_complex c)
{
  /* -i*atanh(i*c) */
  mrb_float x1 = cmath_creal(c);
  mrb_float y1 = cmath_cimag(c);
  mrb_complex c2 = cmath_build_complex(-y1, +x1);
  mrb_complex d2 = cmath_catanh(c2);
  mrb_float x2 = cmath_creal(d2);
  mrb_float y2 = cmath_cimag(d2);
  return cmath_build_complex(+y2, -x2);
}

/* exp(z): return the exponential of z */
DEF_CMATH_METHOD(exp)

/* log(z): return the natural logarithm of z, with branch cut along the negative real axis */
static mrb_value
cmath_log(mrb_state *mrb, mrb_value self) {
  mrb_value z;
  mrb_float base;
  mrb_float real, imag;

  mrb_int n = mrb_get_args(mrb, "o|f", &z, &base);

#ifndef M_E
#define M_E F(exp)(1.0)
#endif

  if (n == 1) base = M_E;
  if (cmath_get_complex(mrb, z, &real, &imag) || real < 0.0) {
    mrb_complex c = cmath_build_complex(real,imag);
    c = cmath_clog(c);
    if (n == 2) c = CXDIVc(c, cmath_clog(cmath_build_complex(base,0)));
    return mrb_complex_new(mrb, cmath_creal(c), cmath_cimag(c));
  }
  if (n == 1) return mrb_float_value(mrb, F(log)(real));
  return mrb_float_value(mrb, F(log)(real)/F(log)(base));
}

/* log10(z): return the base-10 logarithm of z, with branch cut along the negative real axis */
static mrb_value
cmath_log10(mrb_state *mrb, mrb_value self) {
  mrb_value z = mrb_get_arg1(mrb);
  mrb_float real, imag;
  if (cmath_get_complex(mrb, z, &real, &imag) || real < 0.0) {
    mrb_complex c = cmath_build_complex(real,imag);
    c = CXDIVf(cmath_clog(c),log(10));
    return mrb_complex_new(mrb, cmath_creal(c), cmath_cimag(c));
  }
  return mrb_float_value(mrb, F(log10)(real));
}

/* log2(z): return the base-2 logarithm of z, with branch cut along the negative real axis */
static mrb_value
cmath_log2(mrb_state *mrb, mrb_value self) {
  mrb_value z = mrb_get_arg1(mrb);
  mrb_float real, imag;
  if (cmath_get_complex(mrb, z, &real, &imag) || real < 0.0) {
    mrb_complex c = cmath_build_complex(real,imag);
    c = CXDIVf(cmath_clog(c),log(2.0));
    return mrb_complex_new(mrb, cmath_creal(c), cmath_cimag(c));
  }
  return mrb_float_value(mrb, F(log2)(real));
}

/* sqrt(z): return square root of z */
static mrb_value
cmath_sqrt(mrb_state *mrb, mrb_value self) {
  mrb_value z = mrb_get_arg1(mrb);
  mrb_float real, imag;
  if (cmath_get_complex(mrb, z, &real, &imag) || real < 0.0) {
    mrb_complex c = cmath_build_complex(real,imag);
    c = cmath_csqrt(c);
    return mrb_complex_new(mrb, cmath_creal(c), cmath_cimag(c));
  }
  return mrb_float_value(mrb, F(sqrt)(real));
}

/* sin(z): sine function */
DEF_CMATH_METHOD(sin)
/* cos(z): cosine function */
DEF_CMATH_METHOD(cos)
/* tan(z): tangent function */
DEF_CMATH_METHOD(tan)
/* asin(z): arc sine function */
DEF_CMATH_METHOD(asin)
/* acos(z): arc cosine function */
DEF_CMATH_METHOD(acos)
/* atan(z): arg tangent function */
DEF_CMATH_METHOD(atan)
/* sinh(z): hyperbolic sine function */
DEF_CMATH_METHOD(sinh)
/* cosh(z): hyperbolic cosine function */
DEF_CMATH_METHOD(cosh)
/* tanh(z): hyperbolic tangent function */
DEF_CMATH_METHOD(tanh)
/* asinh(z): inverse hyperbolic sine function */
DEF_CMATH_METHOD(asinh)
/* acosh(z): inverse hyperbolic cosine function */
DEF_CMATH_METHOD(acosh)
/* atanh(z): inverse hyperbolic tangent function */
DEF_CMATH_METHOD(atanh)

/* ------------------------------------------------------------------------*/

void
mrb_mruby_cmath_alt_gem_init(mrb_state* mrb)
{
  struct RClass *cmath;
  cmath = mrb_define_module(mrb, "CMath");

  mrb_include_module(mrb, cmath, mrb_module_get(mrb, "Math"));

  mrb_define_module_function(mrb, cmath, "sin", cmath_sin, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "cos", cmath_cos, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "tan", cmath_tan, MRB_ARGS_REQ(1));

  mrb_define_module_function(mrb, cmath, "asin", cmath_asin, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "acos", cmath_acos, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "atan", cmath_atan, MRB_ARGS_REQ(1));

  mrb_define_module_function(mrb, cmath, "sinh", cmath_sinh, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "cosh", cmath_cosh, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "tanh", cmath_tanh, MRB_ARGS_REQ(1));

  mrb_define_module_function(mrb, cmath, "asinh", cmath_asinh, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "acosh", cmath_acosh, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "atanh", cmath_atanh, MRB_ARGS_REQ(1));

  mrb_define_module_function(mrb, cmath, "exp", cmath_exp, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "log", cmath_log, MRB_ARGS_REQ(1)|MRB_ARGS_OPT(1));
  mrb_define_module_function(mrb, cmath, "log2", cmath_log2, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "log10", cmath_log10, MRB_ARGS_REQ(1));
  mrb_define_module_function(mrb, cmath, "sqrt", cmath_sqrt, MRB_ARGS_REQ(1));
}

void
mrb_mruby_cmath_alt_gem_final(mrb_state* mrb)
{
}
