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

    // fft computation (for each block)
    fft(input, spectrum, NFFT);

    // spectral modification (for each block)
    for (i=0; i<NFFT/2+1; i++)
        spectrum[i] = complex_mul(spectrum[i], transfer[i]);

    // ifft computation (for each block)
    ifft(spectrum, output, NFFT);
*/

#include "complex_float32.h"
#include "complex_float64.h"

void fft(float *input, complex_float32 *spectrum, int nfft);
void ifft(complex_float32 *spectrum, float *output, int nfft);
void magnitude(complex_float32 *input, float *result, int n);
void magnitude_db(complex_float32 *input, float *result, int n);
void phase_rad(complex_float32 *input, float *result, int n);

void set_twiddle_table(int max_nfft);
void fft_double(double *input, complex_float64 *spectrum, int n);
void ifft_double(complex_float64 *spectrum, double *output, int n);
void magnitude_double(complex_float64 *input, double *result, int n);
void magnitude_db_double(complex_float64 *input, double *result, int n);
void phase_rad_double(complex_float64 *input, double *result, int n);
int ilog2(int iarg);

#endif
