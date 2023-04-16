#ifndef MicroMath_h
#define MicroMath_h
#include <Arduino.h>

class Complex {
public:
  double re;
  double im;
  Complex(double re, double im);
  Complex(double re);
  Complex(const Complex &a);
  Complex();
  // Addition
  friend Complex operator+(Complex a, Complex b);
  friend Complex operator+(Complex a, double b);
  friend Complex operator+(double a, Complex b);

  Complex &operator+=(const Complex &a);
  Complex &operator+=(const double &b);
  // Subtraction and Negative
  friend Complex operator-(Complex a, Complex b);
  friend Complex operator-(Complex a, double b);
  friend Complex operator-(double a, Complex b);
  friend Complex operator-(Complex a);

  Complex &operator-=(const Complex &a);
  Complex &operator-=(const double &b);
  // Multiplication
  friend Complex operator*(Complex a, Complex b);
  friend Complex operator*(Complex a, double b);
  friend Complex operator*(double a, Complex b);
  // Divison
  friend Complex operator/(Complex a, Complex b);
  friend Complex operator/(Complex a, double b);
  friend Complex operator/(double a, Complex b);
};

class MicroMath {
public:
  Complex I;
  double lambda;
  MicroMath();
  // Complex exponentiation
  Complex c_pow(Complex a, Complex b);
  Complex c_pow(Complex a, double b);
  Complex c_pow(double a, Complex b);
  Complex c_pow(double a, double b);

  Complex c_exp(Complex a);
  Complex c_exp(double a);
  // Sqrt
  Complex c_sqrt(Complex a);
  Complex c_sqrt(double a);
  // Absolute
  double c_abs(Complex a);
  double c_abs(double a);
  // Argument
  double c_arg(Complex a);
  double c_arg(double a);
  // Complex Logarithm
  Complex c_log(Complex a);
  Complex c_log(double a);

  //----------------Trigonometery----------------
  // Sine
  Complex c_sin(Complex a);
  Complex c_sin(double a);
  Complex c_arcsin(Complex a);
  Complex c_arcsin(double a);
  // Cosine
  Complex c_cos(Complex a);
  Complex c_cos(double a);
  Complex c_arccos(Complex a);
  Complex c_arccos(double a);
  // Tangent
  Complex c_tan(Complex a);
  Complex c_tan(double a);
  Complex c_arctan(Complex a);
  Complex c_arctan(double a);
  // Secant
  Complex c_sec(Complex a);
  Complex c_sec(double a);
  Complex c_arcsec(Complex a);
  Complex c_arcsec(double a);
  // Cosecant
  Complex c_csc(Complex a);
  Complex c_csc(double a);
  Complex c_arccsc(Complex a);
  Complex c_arccsc(double a);
  // Cotangent
  Complex c_cot(Complex a);
  Complex c_cot(double a);
  Complex c_arccot(Complex a);
  Complex c_arccot(double a);
  //----------------Hyperbolic Trigonometery----------------
  // Hyperbolic Sine
  Complex c_sinh(Complex a);
  Complex c_sinh(double a);
  Complex c_arcsinh(Complex a);
  Complex c_arcsinh(double a);

  // Hyperbolic Cosine
  Complex c_cosh(Complex a);
  Complex c_cosh(double a);
  Complex c_arccosh(Complex a);
  Complex c_arccosh(double a);

  // Hyperbolic Tangent
  Complex c_tanh(Complex a);
  Complex c_tanh(double a);
  Complex c_arctanh(Complex a);
  Complex c_arctanh(double a);

  // Hyperbolic Cotangent
  Complex c_coth(Complex a);
  Complex c_coth(double a);
  Complex c_arccoth(Complex a);
  Complex c_arccoth(double a);

  // Hyperbolic Secant
  Complex c_sech(Complex a);
  Complex c_sech(double a);
  Complex c_arcsech(Complex a);
  Complex c_arcsech(double a);

  // Hyperbolic Cosecant
  Complex c_csch(Complex a);
  Complex c_csch(double a);
  Complex c_arccsch(Complex a);
  Complex c_arccsch(double a);

  //----------------Special Functions----------------

  // Factorial (Basic)
  Complex factorial(Complex n);
  double factorial(double n);
  int factorial(int n);
  // nCr
  Complex nCr(Complex n, Complex r);
  Complex nCr(Complex n, double r);
  Complex nCr(double n, Complex r);
  double nCr(double n, double r);
  int nCr(int n, int r);
  // Gamma Function
  Complex gamma(Complex z);
  double gamma(double z);
  Complex logGamma(Complex z);
  double logGamma(double z);

  // Chebyshev Polynomials of the First Kind
  Complex chebyshevT(Complex n, Complex x);
  Complex chebyshevT(double n, Complex x);
  Complex chebyshevT(Complex n, double x);
  double chebyshevT(double n, double x);

  // Lambert W
  Complex lambertW(Complex z);
  Complex lambertW(double z);

  // Error Function
  Complex erf(Complex z);
  double erf(double z);

  // Riemann Zeta
  Complex zeta(Complex s);
  double zeta(double s);

  // Log integral
  Complex li(Complex x);
  Complex li(double x);

  // Exponential integral
  Complex ei(Complex x);
  Complex ei(double x);
};

#endif
