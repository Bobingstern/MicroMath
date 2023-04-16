#include "MicroMath.h"

Complex::Complex(double re, double im) {
  this->re = re;
  this->im = im;
}
Complex::Complex(double re) {
  this->re = re;
  this->im = 0;
}
Complex::Complex(const Complex &a) {
  re = a.re;
  im = a.im;
}
Complex::Complex() {
  this->re = 0;
  this->im = 0;
}
// Addition
MicroMath::MicroMath() {
  this->I = Complex(0, 1);
  this->lambda = 0.5772156;
}
Complex operator+(Complex a, Complex b) {
  return Complex(a.re + b.re, a.im + b.im);
}
Complex operator+(Complex a, double b) { return Complex(a.re + b, a.im); }
Complex operator+(double b, Complex a) { return Complex(a.re + b, a.im); }

Complex &Complex::operator+=(const Complex &a) {
  re += a.re;
  im += a.im;
  return *this;
}
Complex &Complex::operator+=(const double &a) {
  re += a;
  return *this;
}
// Subtraction
Complex operator-(Complex a) { return Complex(-a.re, -a.im); }
Complex operator-(Complex a, Complex b) { return a + (-b); }
Complex operator-(Complex a, double b) { return Complex(a.re - b, a.im); }
Complex operator-(double b, Complex a) { return Complex(b - a.re, -a.im); }

Complex &Complex::operator-=(const Complex &a) {
  re -= a.re;
  im -= a.im;
  return *this;
}
Complex &Complex::operator-=(const double &a) {
  re -= a;
  return *this;
}

// Multiplication
Complex operator*(Complex a, Complex b) {
  return Complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
Complex operator*(Complex a, double b) { return Complex(a.re * b, a.im * b); }
Complex operator*(double b, Complex a) { return Complex(a.re * b, a.im * b); }

Complex operator/(Complex a, Complex b) {
  return a * Complex(b.re / (b.re * b.re + b.im * b.im),
                     -b.im / (b.re * b.re + b.im * b.im));
}
Complex operator/(Complex a, double b) { return Complex(a.re / b, a.im / b); }
Complex operator/(double b, Complex a) {
  return Complex(b, 0) / Complex(a.re, a.im);
}

// Exponentiation
Complex MicroMath::c_pow(Complex a_, Complex b_) {
  // double r = sqrt(a.re * a.re + a.im * a.im);
  // double theta = atan2(a.im, a.re);
  // double m = exp(b.re * log(r) - b.im * theta);
  // double theta2 = log(r) * b.re + theta * b.re;
  // return Complex(m * cos(theta2), m * sin(theta2));
  double a = a_.re;
  double b = a_.im;
  double c = b_.re;
  double d = b_.im;

  double f = pow((a * a + b * b), (c / 2)) * exp(-d * atan2(b, a));
  double g = cos(c * atan2(b, a) + d * log(a * a + b * b) / 2);
  double h = sin(c * atan2(b, a) + d * log(a * a + b * b) / 2);
  return f * Complex(g, h);
}
Complex MicroMath::c_pow(Complex a, double b) {
  return c_pow(a, Complex(b, 0));
}
Complex MicroMath::c_pow(double a, Complex b) {
  return c_pow(Complex(a, 0), b);
}
Complex MicroMath::c_pow(double a, double b) {
  return c_pow(Complex(a, 0), Complex(b, 0));
}
// Sqrt
Complex MicroMath::c_sqrt(Complex a) { return c_pow(a, 0.5); }
Complex MicroMath::c_sqrt(double a) { return Complex(sqrt(a), 0); }
// Absolute
double MicroMath::c_abs(Complex a) { return sqrt(a.re * a.re + a.im * a.im); }
double MicroMath::c_abs(double a) { return abs(a); }

// Complex Argument
double MicroMath::c_arg(Complex a) { return atan2(a.im, a.re); }
double MicroMath::c_arg(double a) { return 0; }

// Logarithm
Complex MicroMath::c_log(Complex a) { return Complex(log(c_abs(a)), c_arg(a)); }
Complex MicroMath::c_log(double a) { return c_log(Complex(a, 0)); }
// Exponential
Complex MicroMath::c_exp(Complex a) {
  return Complex(exp(a.re) * cos(a.im), exp(a.re) * sin(a.im));
}
Complex MicroMath::c_exp(double a) { return Complex(exp(a), 0); }
//----------------Trigonometery----------------
// Sine
Complex MicroMath::c_sin(Complex a) {
  return Complex(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}
Complex MicroMath::c_sin(double a) { return Complex(sin(a), 0); }
Complex MicroMath::c_arcsin(Complex a) {
  if (a.re == 1 && a.im == 0)
    return M_PI / 2;
  if (a.re == -1 && a.im == 0)
    return -M_PI / 2;
  return 1 / I * c_log(I * a + c_pow(1 - a * a, 0.5));
}
Complex MicroMath::c_arcsin(double a) { return c_arcsin(Complex(a, 0)); }

// Cosine
Complex MicroMath::c_cos(Complex a) {
  return Complex(cos(a.re) * cosh(a.im), -sin(a.re) * sinh(a.im));
}
Complex MicroMath::c_cos(double a) { return Complex(cos(a), 0); }
Complex MicroMath::c_arccos(Complex a) {
  if (a.re == 1 && a.im == 0)
    return 0;
  if (a.re == -1 && a.im == 0)
    return M_PI;
  return 1 / I * c_log(a + c_pow(a * a - 1, 0.5));
}
Complex MicroMath::c_arccos(double a) { return c_arccos(Complex(a, 0)); }

// Tangent
Complex MicroMath::c_tan(Complex a) {
  return (sin(2 * a.re) + I * sinh(2 * a.im)) /
         (cos(2 * a.re) + cosh(2 * a.im));
}
Complex MicroMath::c_tan(double a) { return Complex(tan(a), 0); }
Complex MicroMath::c_arctan(Complex a) {
  return 1 / (2 * I) * c_log((I - a) / (I + a));
}
Complex MicroMath::c_arctan(double a) { return Complex(atan(a), 0); }

// Secant
Complex MicroMath::c_sec(Complex a) { return 1 / c_cos(a); }
Complex MicroMath::c_sec(double a) { return 1 / cos(a); }
Complex MicroMath::c_arcsec(Complex a) {
  if (a.re == 1 && a.im == 0)
    return 0;
  if (a.re == -1 && a.im == 0)
    return M_PI;
  return 0.5 * M_PI + I * c_log(c_sqrt(1 - 1 / (a * a)) + I / a);
}
Complex MicroMath::c_arcsec(double a) { return Complex(acos(1 / a), 0); }

// Secant
Complex MicroMath::c_csc(Complex a) { return 1 / c_sin(a); }
Complex MicroMath::c_csc(double a) { return 1 / sin(a); }
Complex MicroMath::c_arccsc(Complex a) {
  if (a.re == 1 && a.im == 0)
    return M_PI / 2;
  if (a.re == -1 && a.im == 0)
    return -M_PI / 2;
  return -I * c_log(c_sqrt(1 - 1 / (a * a)) + I / a);
}
Complex MicroMath::c_arccsc(double a) { return Complex(asin(1 / a), 0); }

// Cotangent
Complex MicroMath::c_cot(Complex a) { return 1 / c_tan(a); }
Complex MicroMath::c_cot(double a) { return 1 / tan(a); }
Complex MicroMath::c_arccot(Complex a) { return c_arctan(1 / a); }
Complex MicroMath::c_arccot(double a) { return Complex(atan(1 / a), 0); }

//----------------Hyperbolic Trigonometery----------------
// Hyperbolic Sine
Complex MicroMath::c_sinh(Complex a) {
  return Complex(sinh(a.re) * cos(a.im), cosh(a.re) * sin(a.im));
}
Complex MicroMath::c_sinh(double a) { return Complex(sinh(a), 0); }
Complex MicroMath::c_arcsinh(Complex a) { return c_log(a + c_sqrt(1 + a * a)); }
Complex MicroMath::c_arcsinh(double a) { return Complex(log(a + sqrt(1 + a * a)), 0); }

// Hyperbolic Cosine
Complex MicroMath::c_cosh(Complex a) {
  return Complex(cosh(a.re) * cos(a.im), sinh(a.re) * sin(a.im));
}
Complex MicroMath::c_cosh(double a) { return Complex(cosh(a), 0); }
Complex MicroMath::c_arccosh(Complex a) {
  return c_log(a + c_sqrt(a - 1) * c_sqrt(1 + a));
}
Complex MicroMath::c_arccosh(double a) { return c_arccosh(Complex(a, 0)); }

// Hyperbolic Tangent
Complex MicroMath::c_tanh(Complex a) {
  return (sinh(2 * a.re) + I * sin(2 * a.im)) /
         (cosh(2 * a.re) + cos(2 * a.im));
}
Complex MicroMath::c_tanh(double a) { return Complex(tanh(a), 0); }
Complex MicroMath::c_arctanh(Complex a) {
  return 0.5 * (c_log(1 + a) - c_log(1 - a));
}
Complex MicroMath::c_arctanh(double a) { return c_arctanh(Complex(a, 0)); }

// Hyperbolic Cotangent
Complex MicroMath::c_coth(Complex a) { return 1 / c_tanh(a); }
Complex MicroMath::c_coth(double a) { return 1 / c_tanh(a); }
Complex MicroMath::c_arccoth(Complex a) { return c_arctanh(1 / a); }
Complex MicroMath::c_arccoth(double a) { return c_arctanh(1 / a); }

// Hyperbolic Secant
Complex MicroMath::c_sech(Complex a) { return 1 / c_cosh(a); }
Complex MicroMath::c_sech(double a) { return 1 / c_cosh(a); }
Complex MicroMath::c_arcsech(Complex a) { return c_arccosh(1 / a); }
Complex MicroMath::c_arcsech(double a) { return c_arccosh(1 / a); }

// Hyperbolic Cosecant
Complex MicroMath::c_csch(Complex a) { return 1 / c_sinh(a); }
Complex MicroMath::c_csch(double a) { return 1 / c_sinh(a); }
Complex MicroMath::c_arccsch(Complex a) { return c_arcsinh(1 / a); }
Complex MicroMath::c_arccsch(double a) { return c_arcsinh(1 / a); }

//----------------Special Functions----------------

// Factorial
Complex MicroMath::factorial(Complex n) { return gamma(n + 1); }
double MicroMath::factorial(double n) { return gamma(n + 1); }
int MicroMath::factorial(int n) {
  int y = 1;
  for (int i = 2; i <= n; i++) {
    y *= i;
  }
  return y;
}
// nCr
Complex MicroMath::nCr(Complex n, Complex r) {
  return gamma(n + 1) / (gamma(r + 1) * gamma(n - r + 1));
}
Complex MicroMath::nCr(Complex n, double r) {
  return gamma(n + 1) / (gamma(r + 1) * gamma(n - r + 1));
}
Complex MicroMath::nCr(double n, Complex r) {
  return gamma(n + 1) / (gamma(r + 1) * gamma(n - r + 1));
}
double MicroMath::nCr(double n, double r) {
  return gamma(n + 1) / (gamma(r + 1) * gamma(n - r + 1));
}
int MicroMath::nCr(int n, int r) {
  long long y = 1;
  for (int i = 1; i <= r; i++) {
    y *= n - i + 1;
    y /= i;
  }
  return y;
}
// Gamma
Complex MicroMath::gamma(Complex z) {
  Complex y;
  if (z.re < 0.5) {
    y = M_PI / (c_sin(M_PI * z) * gamma(1 - z));
  } else {
    y = c_sqrt(2 * M_PI / z) *
        c_pow(1 / M_E * (z + 1 / (12 * z - 1 / (10 * z))), z);
  }
  return y;
}
double MicroMath::gamma(double z) {
  double y;
  if (z < 0.5) {
    y = M_PI / (sin(M_PI * z) * gamma(1 - z));
  } else {
    y = sqrt(2 * M_PI / z) *
        pow(1 / M_E * (z + 1 / (12 * z - 1 / (10 * z))), z);
  }
  return y;
}
Complex MicroMath::logGamma(Complex z) {
  Complex y;
  if (z.re < 0.5) {
    y = c_log(M_PI / (c_sin(M_PI * z))) - logGamma(1 - z);
  } else {
    y = 0.5 * (log(2 * M_PI) - c_log(z)) +
        z * (c_log(z + 1 / (12 * z - 1 / (10 * z))) - 1);
  }
  return y;
}
double MicroMath::logGamma(double z) {
  double y;
  if (z < 0.5) {
    y = log(M_PI / (sin(M_PI * z))) - logGamma(1 - z);
  } else {
    y = 0.5 * (log(2 * M_PI) - log(z)) +
        z * (log(z + 1 / (12 * z - 1 / (10 * z))) - 1);
  }
  return y;
}

// Chebyshev Polynomials of the First Kind
Complex MicroMath::chebyshevT(Complex n, Complex x) {
  if (c_abs(x) <= 1)
    return c_cos(n * c_arccos(x));
  return 0.5 *
         (c_pow(x - c_sqrt(x * x - 1), n) + c_pow(x + c_sqrt(x * x - 1), n));
}
Complex MicroMath::chebyshevT(double n, Complex x) {
  if (c_abs(x) <= 1)
    return c_cos(n * c_arccos(x));
  return 0.5 *
         (c_pow(x - c_sqrt(x * x - 1), n) + c_pow(x + c_sqrt(x * x - 1), n));
}
Complex MicroMath::chebyshevT(Complex n, double x) {
  if (c_abs(x) <= 1)
    return c_cos(n * c_arccos(x));
  return 0.5 *
         (c_pow(x - c_sqrt(x * x - 1), n) + c_pow(x + c_sqrt(x * x - 1), n));
}
double MicroMath::chebyshevT(double n, double x) {
  if (c_abs(x) <= 1)
    return (c_cos(n * c_arccos(x))).re;
  return (0.5 *
          (c_pow(x - c_sqrt(x * x - 1), n) + c_pow(x + c_sqrt(x * x - 1), n)))
      .re;
}

// Lambert W
Complex MicroMath::lambertW(Complex z) {
  Complex cq = c_sqrt(1 + M_E * z);
  Complex w = M_E * z / (1 + M_E * z + cq) * c_log(1 + cq);
  Complex old(w);
  // Max 10 Halley iterations
  for (int i = 0; i < 10; i++) {
    Complex ew = c_exp(w);
    w = w -
        (w * ew - z) / (ew * (w + 1) - ((w + 2) * (w * ew - z)) / (2 * w + 2));
    if (abs(c_abs(w) - c_abs(old)) < 1e-6) {
      break;
    }
    old = w;
  }
  return w;
}
Complex MicroMath::lambertW(double z) {
  double cq = sqrt(1 + M_E * z);
  double w = M_E * z / (1 + M_E * z + cq) * log(1 + cq);
  double old = w;
  // Max 10 Halley iterations
  for (int i = 0; i < 10; i++) {
    double ew = exp(w);
    w = w -
        (w * ew - z) / (ew * (w + 1) - ((w + 2) * (w * ew - z)) / (2 * w + 2));
    if (abs(w - old) < 1e-6) {
      break;
    }
    old = w;
  }
  return w;
}

// Error Function
Complex MicroMath::erf(Complex z) {
  Complex a =
      erf(z.re) + exp(-z.re * z.re) / (2 * M_PI * z.re) *
                      (1 - cos(2 * z.re * z.im) + I * sin(2 * z.re * z.im));
  Complex s = 0;
  Complex old = 0;
  for (int k = 1; k < 16; k++) {
    Complex f = 2 * z.re * (1 - cos(2 * z.re * z.im) * cosh(k * z.im)) +
                k * sin(2 * z.re * z.im) * sinh(k * z.im);
    Complex g = 2 * z.re * sin(2 * z.re * z.im) * cosh(k * z.im) +
                k * cos(2 * z.re * z.im) * sinh(k * z.im);

    s += exp(-k * k * 0.25) / (k * k + 4 * z.re * z.re) * (f + I * g);
    if (abs(c_abs(s) - c_abs(old)) < 1e-6) {
      break;
    }
    old = s;
  }
  return a + 2 * exp(-z.re * z.re) / M_PI * s;
}
double MicroMath::erf(double z) {
  if (z <= 0) {
    return -erf(-z);
  }
  double t = 1 / (1 + 0.3275911 * z);

  return 1 - exp(-z * z) * (0.254829592 * t - 0.284496736 * pow(t, 2) +
                            1.421413741 * pow(t, 3) - 1.453152027 * pow(t, 4) +
                            1.061405429 * pow(t, 5));
}

// Riemann Zeta
Complex MicroMath::zeta(Complex s) {
  Complex z;
  Complex y;
  if (s.re < 1) {
    return zeta(1 - s) * 2 * gamma(1 - s) / c_pow(2 * M_PI, 1 - s) *
           c_cos(M_PI * (1 - s) / 2);
  } else {
    int n = floor(1.3 * 6 + 0.9 * abs(s.im));
    if (n > 16)
      n = 16;
    for (int k = 1; k <= n; k++) {
      double dk = 0;
      for (int j = k; j <= n; j++) {
        dk += (factorial(n + j - 1) * pow(4, j)) /
              (factorial(n - j) * factorial(2 * j));
      }
      dk *= n;
      z += pow(-1, k - 1) * dk / c_pow(k, s);
    }
    double dk0 = 0;
    for (int j = 0; j <= n; j++) {
      dk0 += (factorial(n + j - 1) * pow(4, j)) /
             (factorial(n - j) * factorial(2 * j));
    }
    dk0 *= n;
    y = 1 / (dk0 * (1 - c_pow(2, 1 - s))) * z;
  }
  return y;
}
double MicroMath::zeta(double s) { return zeta(Complex(s, 0)).re; }

// Log Integral
Complex MicroMath::li(Complex x) {
  Complex s;
  Complex old;
  for (int n = 1; n < 16; n++) {
    double nest = 0;
    for (int k = 0; k <= floor((n - 1) * 0.5); k++) {
      nest += 1 / ((double)(2 * k + 1));
    }
    s += nest * pow(-1, n - 1) * c_pow(c_log(x), n) /
         (double)(factorial(n) * pow(2, n - 1));

    if (abs(c_abs(s) - c_abs(old)) < 1e-6) {
      break;
    }
    old = s;
  }
  return lambda + c_log(c_log(x)) + c_sqrt(x) * s;
}
Complex MicroMath::li(double x) { return li(Complex(x, 0)); }

// Exponential Integral
Complex MicroMath::ei(Complex x) { return li(c_exp(x)); }
Complex MicroMath::ei(double x) { return li(exp(x)); }

// Lanczos approximation. Overkill for microcontrollers
// static const int MicroMath_GammaPlen = 8;
// double MicroMath_GammaP[MicroMath_GammaPlen] = {
//     676.5203681218851,     -1259.1392167224028,  771.32342877765313,
//     -176.61502916214059,   12.507343278686905,   -0.13857109526572012,
//     9.9843695780195716e-6, 1.5056327351493116e-7};
//  Complex MicroMath::gamma(Complex z) {
//    Complex y;
//    if (z.re < 0.5) {
//      y = M_PI / (c_sin(M_PI * z) * gamma(1 - z));
//    } else {
//      Complex x(1, 0);
//      for (int i = 0; i < MicroMath_GammaPlen; i++) {
//        x += MicroMath_GammaP[i] / (z - 1 + i + 1);
//      }
//      Complex t = z - 1 + 8 - 0.5;
//      y = sqrt(2 * M_PI) * c_pow(t, z - 1 + 0.5) * c_exp(-t) * x;
//    }
//    return y;
//  }
//  double MicroMath::gamma(double z) {
//    double y;
//    if (z < 0.5) {
//      y = M_PI / (sin(M_PI * z) * gamma(1 - z));
//    } else {
//      double x = 1;
//      for (int i = 0; i < MicroMath_GammaPlen; i++) {
//        x += MicroMath_GammaP[i] / (z - 1 + i + 1);
//      }
//      double t = z - 1 + 8 - 0.5;
//      y = sqrt(2 * M_PI) * pow(t, z - 1 + 0.5) * exp(-t) * x;
//    }
//    return y;
//  }