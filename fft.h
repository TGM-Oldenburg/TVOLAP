/*----------------------------------------------------------------------------*\
| Header of fft.cpp, for explanation and exampple, see cpp-file.               |
|                                                                              |
| Author: (c) Uwe Simmer                                  June 1988 - Nov 2012 |
| LGPL Release: May 2017, License see end of file                              |
\*----------------------------------------------------------------------------*/

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

/*------------------------------License---------------------------------------*\
| Copyright (c) 1988-2012 Uwe Simmer                                           |
|                                                                              |
| This program is free software: you can redistribute it and/or modify         |
| it under the terms of the GNU Lesser General Public License as published by  |
| the Free Software Foundation, either version 3 of the License, or            |
| (at your option) any later version.                                          |
|                                                                              |
| This program is distributed in the hope that it will be useful,              |
| but WITHOUT ANY WARRANTY; without even the implied warranty of               |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                |
| GNU Lesser General Public License for more details.                          |
|                                                                              |
| You should have received a copy of the GNU Lesser General Public License     |
| along with this program. If not, see <http://www.gnu.org/licenses/>.         |
\*----------------------------------------------------------------------------*/
