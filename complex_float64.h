#ifndef _COMPLEX_FLOAT64
#define _COMPLEX_FLOAT64

typedef struct {
    double re,im;
} complex_float64;

#include <math.h>

//------------------------------------------------------------------------------

inline complex_float64 complex(double real_part, double imag_part)
{
    complex_float64 z;

    z.re = real_part;
    z.im = imag_part;

    return z;
}

//------------------------------------------------------------------------------

inline double complex_real(complex_float64 x)
{
    return x.re;
}

//------------------------------------------------------------------------------

inline double complex_imag(complex_float64 x)
{
    return x.im;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_add(complex_float64 x, complex_float64 y)
{
    complex_float64 z;

    z.re = x.re + y.re;
    z.im = x.im + y.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_sub(complex_float64 x, complex_float64 y)
{
    complex_float64 z;

    z.re = x.re - y.re;
    z.im = x.im - y.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_neg(complex_float64 x)
{
    complex_float64 z;

    z.re = -x.re;
    z.im = -x.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_conj(complex_float64 x)
{
    complex_float64 z;

    z.re =  x.re;
    z.im = -x.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_mul(complex_float64 x, complex_float64 y)
{
    complex_float64 z;

    z.re = (x.re * y.re)
         - (x.im * y.im);

    z.im = (x.re * y.im)
         + (x.im * y.re);

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_mul(complex_float64 x, double y)
{
    complex_float64 z;

    z.re = x.re * y;
    z.im = x.im * y;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_div(complex_float64 numerator,
                                  complex_float64 denominator)
{
    double divisor;
    complex_float64 quotient;

    divisor = denominator.re * denominator.re
            + denominator.im * denominator.im;

    quotient.re = (numerator.re * denominator.re
                +  numerator.im * denominator.im) / divisor;

    quotient.im = (numerator.im * denominator.re
                -  numerator.re * denominator.im) / divisor;

    return quotient;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_div(complex_float64 numerator, double denominator)
{
    complex_float64 quotient;

    quotient.re = numerator.re / denominator;
    quotient.im = numerator.im / denominator;

    return quotient;
}

//------------------------------------------------------------------------------

inline complex_float64 complex_inv(complex_float64 x)
{
    double divisor;
    complex_float64 inverse;

    divisor = x.re * x.re + x.im * x.im;

    inverse.re =  x.re / divisor;
    inverse.im = -x.im / divisor;

    return inverse;
}

//------------------------------------------------------------------------------

inline double complex_abs(complex_float64 x)
{
    return sqrt(x.re * x.re + x.im * x.im);
}

//------------------------------------------------------------------------------

inline double complex_abs_squared(complex_float64 x)
{
    return (x.re * x.re + x.im * x.im);
}

//------------------------------------------------------------------------------

inline double complex_angle(complex_float64 x)
{
    return atan2(x.im, x.re);
}

//------------------------------------------------------------------------------

inline complex_float64 complex_sqrt(complex_float64 x)
{
    double x_abs, temp;
    complex_float64 z;

    x_abs = complex_abs(x);
    temp = (x_abs + x.re) / 2.;
    if (temp < 0.)
        temp = 0.;

    z.re = sqrt(temp);

    temp = (x_abs - x.re) / 2.;
    if (temp < 0.)
        temp = 0.;

    if (x.im > 0.)
        z.im =  sqrt(temp);
    else
        z.im = -sqrt(temp);

    return z;
}

//--------------------- License ------------------------------------------------

// Copyright (c) 2012 Uwe Simmer

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, 
// including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, 
// and to permit persons to whom the Software is furnished to do so, 
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//------------------------------------------------------------------------------


#endif  // _COMPLEX_FLOAT64
