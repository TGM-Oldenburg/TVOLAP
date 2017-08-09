/*----------------------------------------------------------------------------*\
| Header of fft.cpp, for explanation and exampple, see cpp-file.               |
|                                                                              |
| Author: (c) Uwe Simmer                                  June 1988 - Nov 2012 |
| MIT Release: Aug 2017, License see end of file                               |
\*----------------------------------------------------------------------------*/

#ifndef NULL
#define NULL 0x0
#endif

#ifndef _FFT
#define _FFT

#include "complex_float32.h"
#include "complex_float64.h"

void rfft(float *input, complex_float32 *spectrum, int nfft);
void irfft(complex_float32 *spectrum, float *output, int nfft);
void cfft(complex_float32 *x, int nfft);
void icfft(complex_float32 *x, int nfft);

void magnitude(complex_float32 *input, float *result, int n);
void magnitude_db(complex_float32 *input, float *result, int n);
void phase_rad(complex_float32 *input, float *result, int n);

void set_twiddle_table(int max_nfft);
void rfft_double(double *input, complex_float64 *spectrum, int n);
void irfft_double(complex_float64 *spectrum, double *output, int n);
void cfft_double(complex_float64 *x, int nfft);
void icfft_double(complex_float64 *x, int nfft);

void magnitude_double(complex_float64 *input, double *result, int n);
void magnitude_db_double(complex_float64 *input, double *result, int n);
void phase_rad_double(complex_float64 *input, double *result, int n);
int ilog2(int iarg);

int table_get_nfft();
complex_float64 *table_get_twiddle_factor();
double *table_get_cos_half();

#endif

/*------------------------------License----------------------------------------*\
| Copyright (c) 2012-2017 Hagen Jaeger                           				|
|																				|
| Permission is hereby granted, free of charge, to any person obtaining a 		|
| copy of this software and associated documentation files (the "Software"), 	|
| to deal in the Software without restriction, including without limitation 	|
| the rights to use, copy, modify, merge, publish, distribute, sublicense, 		|
| and/or sell copies of the Software, and to permit persons to whom the 		|
| Software is furnished to do so, subject to the following conditions:			|
|																				|
| The above copyright notice and this permission notice shall be included 		|
| in all copies or substantial portions of the Software.						|
|																				|
| THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 		|
| OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 	|
| FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 		|
| THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 	|
| LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 		|
| FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 			|
| DEALINGS IN THE SOFTWARE.														|
|																				|
| https://opensource.org/licenses/mit-license.php								|
\*-----------------------------------------------------------------------------*/
