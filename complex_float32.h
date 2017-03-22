#ifndef _COMPLEX_FLOAT32
#define _COMPLEX_FLOAT32

typedef struct {
    float re,im;
} complex_float32;

#include <math.h>

//------------------------------------------------------------------------------

inline complex_float32 complex(float real_part, float imag_part)
{
    complex_float32 z;

    z.re = real_part;
    z.im = imag_part;

    return z;
}

//------------------------------------------------------------------------------

inline float complex_real(complex_float32 x)
{
    return x.re;
}

//------------------------------------------------------------------------------

inline float complex_imag(complex_float32 x)
{
    return x.im;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_add(complex_float32 x, complex_float32 y)
{
    complex_float32 z;

    z.re = x.re + y.re;
    z.im = x.im + y.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_sub(complex_float32 x, complex_float32 y)
{
    complex_float32 z;

    z.re = x.re - y.re;
    z.im = x.im - y.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_neg(complex_float32 x)
{
    complex_float32 z;

    z.re = -x.re;
    z.im = -x.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_conj(complex_float32 x)
{
    complex_float32 z;

    z.re =  x.re;
    z.im = -x.im;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_mul(complex_float32 x, complex_float32 y)
{
    complex_float32 z;

    z.re = (x.re * y.re)
         - (x.im * y.im);

    z.im = (x.re * y.im)
         + (x.im * y.re);

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_mul(complex_float32 x, float y)
{
    complex_float32 z;

    z.re = x.re * y;
    z.im = x.im * y;

    return z;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_div(complex_float32 numerator,
                                 complex_float32 denominator)
{
    float divisor;
    complex_float32 quotient;

    divisor = denominator.re * denominator.re
            + denominator.im * denominator.im;

    quotient.re = (numerator.re * denominator.re
                +  numerator.im * denominator.im) / divisor;

    quotient.im = (numerator.im * denominator.re
                -  numerator.re * denominator.im) / divisor;

    return quotient;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_div(complex_float32 numerator, float denominator)
{
    complex_float32 quotient;

    quotient.re = numerator.re / denominator;
    quotient.im = numerator.im / denominator;

    return quotient;
}

//------------------------------------------------------------------------------

inline complex_float32 complex_inv(complex_float32 x)
{
    float divisor;
    complex_float32 inverse;

    divisor = x.re * x.re + x.im * x.im;

    inverse.re =  x.re / divisor;
    inverse.im = -x.im / divisor;

    return inverse;
}

//------------------------------------------------------------------------------

inline float complex_abs(complex_float32 x)
{
    return sqrtf(x.re * x.re + x.im * x.im);
}

//------------------------------------------------------------------------------

inline float complex_abs_squared(complex_float32 x)
{
    return (x.re * x.re + x.im * x.im);
}

//------------------------------------------------------------------------------

inline float complex_angle(complex_float32 x)
{
    return atan2f(x.im, x.re);
}

//------------------------------------------------------------------------------

inline complex_float32 complex_sqrt(complex_float32 x)
{
    float x_abs, temp;
    complex_float32 z;

    x_abs = complex_abs(x);
    temp = (x_abs + x.re) / 2.0f;
    if (temp < 0.0f)
        temp = 0.0f;

    z.re = sqrtf(temp);

    temp = (x_abs - x.re) / 2.0f;
    if (temp < 0.0f)
        temp = 0.0f;

    if (x.im > 0.0f)
        z.im =  sqrtf(temp);
    else
        z.im = -sqrtf(temp);

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

#endif  // _COMPLEX_FLOAT32
