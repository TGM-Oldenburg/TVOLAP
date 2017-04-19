#ifndef _FFT
#define _FFT

/*
In-place radix-2 decimation-in-time FFT

Example:

#include "fft.h"
#define NFFT 16384
#define MAX_FFT 65536

    // declarations
    complex_float32 spectrum[NFFT/2+1];
    complex_float32 transfer[NFFT/2+1];
    float input[NFFT];
    float output[NFFT];

    // initialization (once)
    set_twiddle_table(MAX_FFT);

    // real fft computation (for each block)
    rfft(input, spectrum, NFFT);

    // spectral modification (for each block)
    for (i=0; i<NFFT/2+1; i++)
        spectrum[i] = complex_mul(spectrum[i], transfer[i]);

    // inverse real fft computation (for each block) 
    irfft(spectrum, output, NFFT);
*/

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
