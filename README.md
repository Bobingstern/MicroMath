# MicroMath
Full complex support and a bunch of special math functions for Arduino Microcontrollers!


## What is this?
MicroMath is an Arduino library for... Micro Math! It implements complete complex number support and a variety of special math functions. It's also super simple and intuitive to use!

## Installation
MicroMath is not on the Arduino Library Manager yet (it will be soon!) so you will hav to install it manually. Don't worry though, it's super simple! Just grab the `MicroMath.h` and `MicroMath.cpp` files from the latest release or GitHub repo and drag them into your project folder. Your project structure should look something like this:
```
YourProject
|	YourProject.ino
|	MicroMath.h
|	MicroMath.cpp
```

Now you can include the library with `#include "MicroMath.cpp"` at the top of your program.

## Documentation
Here are the docs for how to use MicroMath and its functionality

### Complex
MicroMath comes with complete complex support that is simple and easy to use. The `Complex` class comes with several overloads to create a complex number:
```cpp
#include "MicroMath.h"

Complex a; //0+0i
Complex b(1); //1+0i
Comlex c(2, 3); //2+3i
Complex d(c); //2+3i (Copy constructor)
```
You can access the real and imaginary components of complex numbers using `a.re` and `a.im`. Complex numbers can be added to other complex numbers or doubles through operator overloading:

```cpp
Complex a(2, 3); //2+3i
Complex b(3, 2); //3+2i
Complex c = a + b; //5+5i
Complex d = a - b; //-1+i
Complex e = a * b; //0+13i
Complex f = a / b; //0.9230... + 0.384... i

a += b; //a = 5+5i
//etc..
```

### Standard functions
MicroMath provides all the standard functions in the `math.h` library and more for all Complex valued inputs. This includes all normal and hyperbolic trigonometry and their inverses, log, exp, pow, sqrt, abs and arg. All the complex functions are prefixed with `c_` to avoid confusion with the built-in ones.

#### Power
```cpp
Complex c_pow(Complex a, Complex b);
Complex c_pow(Complex a, double b);
Complex c_pow(double a, Complex b);
Complex c_pow(double a, double b);
```
Calculates $a^b$ for complex or real inputs.

#### Exp
```cpp
Complex c_exp(Complex x);
Complex c_exp(double x);
```
Calculates $e^x$ for complex or real inputs.

#### Log
```cpp
Complex c_log(Complex x);
Complex c_log(double x);
```
Calculates the natural logarithm $\log(x)$ for complex or real inputs.

#### Square root
```cpp
Complex c_sqrt(Complex x);
Complex c_sqrt(double x);
```
Calculates $\sqrt{x}$ for complex or real inputs.

#### Absolute
```cpp
double c_abs(Complex a);
double c_abs(double a);
```
Calculates $|a|$ for complex or real inputs. The absolute of a complex number is its magnitude $|x+yi| = \sqrt{x^2 + y^2}$

#### Argument
```cpp
double c_arg(Complex x);
double c_arg(double x);
```
Calculates $atan2(y,x)$ of $x+yi$ for complex or real inputs.

### Trigonometry
All the trig functions are included in MicroMath. They are all prefixed with `c_` and use the `arc` prefix to denote inverse. I will not document them all explicitly since they are all self explanatory. Here is a section from the header file containing them all:
```cpp
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
```

### Special Functions
MicroMath contains a handful of common special math functions that are used in a variety of fields. They were implemented with microcontrollers in mind so they are not 100% accurate but are designed to be relatively fast. None of them use precomputed arrays of values since I thought that would be overkill.

#### Factorial
```cpp
Complex factorial(Complex n);
double factorial(double n);
int factorial(int n);
```
Calculates the factorial $n!$ of any real or complex value. Makes use of the gamma function (see below) to extend the definition to the entire complex plane.

#### nCr
```cpp
Complex nCr(Complex n, Complex r);
Complex nCr(Complex n, double r);
Complex nCr(double n, Complex r);
double nCr(double n, double r);
int nCr(int n, int r);
```
Calculates the binomial coefficient $n \choose r$ of any 2 real or complex values. Makes use of the gamma function (see below) to extend the definition to the entire complex plane.

#### Gamma Function
```cpp
Complex gamma(Complex z);
double gamma(double z);
```
Calculates the Gamma function $\Gamma(z)$ for any real or complex number.

#### Log Gamma Function
```cpp
Complex logGamma(Complex z);
double logGamma(double z);
```
Calculates the Natural logarithm of the Gamma function $\log \Gamma(z)$ for any real or complex number.

#### Cheybyshev Polynomials of the First Kind
```cpp
Complex chebyshevT(Complex n, Complex x);
Complex chebyshevT(double n, Complex x);
Complex chebyshevT(Complex n, double x);
double chebyshevT(double n, double x);
```
Calculates Chebyshev polynomials $T_n(x)$ for all real and complex values.

#### Lambert W Function
```cpp
Complex lambertW(Complex z);
Complex lambertW(double z);
```
Calculates the Lambert W function defined to be the inverse of the function $f(x)=xe^x$, meaning that $W(x) = f^{-1}(x)$. Valid for all real and complex values.

#### Gaussian Error Function
```cpp
Complex erf(Complex z);
double erf(double z);
```
Calculates the error function $\frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2}dt$ for all real and complex values. This function is more computationally expensive for complex numbers with imaginary components farther from 0.

#### Riemann Zeta function
```cpp
Complex zeta(Complex s);
double zeta(double s);
```
Calculates the Riemann Zeta function $\zeta(s) = \sum_{n=1}^{\infty} \frac{1}{n^s}$ for all complex and real values making use of Riemann's Functional equation to analytically continue the function to the entire complex plane. This function is computationall expensive for large imaginary components of `s`.

#### Logarithmic Integral Function
```cpp
Complex li(Complex x);
Complex li(double x);
```
Calculates the logarthmic integral $li(x) = \int_0^x \frac{dt}{\log t}$ for all real and complex values.

#### Exponential Integral Function
```cpp
Complex ei(Complex x);
Complex ei(double x);
```
Calculates the exponential integral $Ei(x) = \int_{-x}^{\infty} \frac{e^{-t}}{t}$ for all real and complex values. Makes use of the logarithmic integral.

## Contributing
This library is still in its early stages and is always open to contribution! If you find any bugs or issues with it, dont hesitate to open an issue or pr. I'm currently looking for ways to implement Bessel functions and other functions so you have any ideas, please let me know!
