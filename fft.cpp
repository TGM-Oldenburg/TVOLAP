/*----------------------------------------------------------------------------*\
| In-place radix-2 decimation in time FFT                                      |
|                                                                              |
|   Example:                                                                   |
|                                                                              |
|   #include "fft.h"                                                           |
|   #define NFFT 16384                                                         |
|   #define MAX_FFT 65536                                                      |
|                                                                              |
|   // declarations                                                            |
|   complex_float32 spectrum[NFFT/2+1];                                        |
|   complex_float32 transfer[NFFT/2+1];                                        |
|   float input[NFFT];                                                         |
|   float output[NFFT];                                                        |
|                                                                              |
|   // initialization (once)                                                   |
|   set_twiddle_table(MAX_FFT);                                                |
|                                                                              |
|   // real fft computation (for each block)                                   |
|   rfft(input, spectrum, NFFT);                                               |
|                                                                              |
|   // spectral modification (for each block)                                  |
|   for (i=0; i<NFFT/2+1; i++)                                                 |
|       spectrum[i] = complex_mul(spectrum[i], transfer[i]);                   |
|                                                                              |
|   // inverse real fft computation (for each block)                           |
|   irfft(spectrum, output, NFFT);                                             |
|                                                                              |
| Author: (c) Uwe Simmer                                  June 1988 - Nov 2012 |
| MIT Release: Aug. 2017, License see end of file                              |
\*----------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft.h"

typedef struct {
    int nfft;
    complex_float64 *twiddle_factor;
    double *cos_half;
} fft_table;

static fft_table table;

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#ifndef NULL
#define NULL 0x0
#endif

//------------------------------------------------------------------------------

void set_twiddle_table(int max_nfft)
{
    int i;

    if (ilog2(max_nfft) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", max_nfft);
        return;
    }

    // real FFT by half-length complex FFT
    table.nfft = max_nfft / 2;

    if (table.twiddle_factor)
        free(table.twiddle_factor);
    table.twiddle_factor = (complex_float64 *) malloc(table.nfft/2*sizeof(complex_float64));
    if (table.twiddle_factor == NULL)
    {
        printf("Error: Insufficient memory.\n");
        return;
    }

    // compute exp(-jw) table for complex fft core
    for (i=0; i<table.nfft/2; i++)
    {
        table.twiddle_factor[i].re = +cos(2. * M_PI * i / table.nfft);
        table.twiddle_factor[i].im = -sin(2. * M_PI * i / table.nfft);
    }

    if (table.cos_half)
        free(table.cos_half);
    table.cos_half = (double *) malloc(table.nfft/2*sizeof(double));
    if (table.cos_half == NULL)
    {
        printf("Error: Insufficient memory.\n");
        return;
    }

    // compute cos table for real fft
    for (i=0; i<table.nfft/2; i++)
    {
        table.cos_half[i] = cos(M_PI * i / table.nfft);
    }
}

//------------------------------------------------------------------------------

int table_get_nfft()
{
    return table.nfft;
}

//------------------------------------------------------------------------------

complex_float64 *table_get_twiddle_factor()
{
    return table.twiddle_factor;
}

//------------------------------------------------------------------------------

double *table_get_cos_half()
{
    return table.cos_half;
}

//------------------------------------------------------------------------------

void fft_core(complex_float32 *x, complex_float64 *w, int nstride, int nfft)
{
    int i, ig, j, k, l, ngroups, nbutterflies;

    complex_float32 ctemp;

    //----- Bit Reverse Section ------------------------------------------------

    j = 0;
    for (i=0; i<nfft-1; i++)
    {
        if (i<j)
        {
            ctemp.re = x[j].re;
            ctemp.im = x[j].im;

            x[j].re = x[i].re;
            x[j].im = x[i].im;

            x[i].re = ctemp.re;
            x[i].im = ctemp.im;
        }

        k = nfft / 2;
        while (k <= j)
        {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    //----- First Stage --------------------------------------------------------

    i = 0;
    j = 1;
    ngroups = nfft >> 1;

    for (ig=0; ig<ngroups; ig++)
    {
        ctemp.re = x[j].re;
        ctemp.im = x[j].im;

        x[j].re = x[i].re - ctemp.re;
        x[j].im = x[i].im - ctemp.im;

        x[i].re = x[i].re + ctemp.re;
        x[i].im = x[i].im + ctemp.im;

        i = i + 2;
        j = j + 2;
    }

    //----- Final log2(nfft)-1 Stages ------------------------------------------

    nbutterflies = 2;

    for (ngroups = nfft >> 2; ngroups > 0; ngroups >>= 1)
    {
        i = 0;
        j = nbutterflies;

        for (ig=0; ig<ngroups; ig++)
        {
            k = 0;
            for (l=0; l<nbutterflies; l++)
            {
                // Butterfly
                ctemp.re = (float) (x[j].re * w[k].re - x[j].im * w[k].im);
                ctemp.im = (float) (x[j].im * w[k].re + x[j].re * w[k].im);

                x[j].re = x[i].re - ctemp.re;
                x[j].im = x[i].im - ctemp.im;

                x[i].re = x[i].re + ctemp.re;
                x[i].im = x[i].im + ctemp.im;

                k += (ngroups << nstride);
                i++;
                j++;
            }

            i = i + nbutterflies;
            j = j + nbutterflies;
        }

        nbutterflies <<= 1;
    }
}

//------------------------------------------------------------------------------

void rfft(float *input, complex_float32 *spectrum, int n)
{
    int i, j, k, nfft, nstride;
    float tr, ti, rs, is, rd, id, rp, ip, ci, cj;
    complex_float32 *x;
    complex_float64 *w;
    double *cos2table;

    nfft = n/2;

    if (nfft == 0 || input == NULL || spectrum == NULL)
        return;

    if (ilog2(n) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", n);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(n);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(n);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2*k, 2*nfft);
        return;
    }

    // copy memory if not in-place
    if (input != (float *)spectrum)
    {
        for (i=0; i<n/2; i++)
        {
            spectrum[i].re = input[2*i + 0];
            spectrum[i].im = input[2*i + 1];
        }
    }

    x = spectrum;
    w = table.twiddle_factor;
    cos2table = table.cos_half;

    fft_core(x, w, nstride, nfft);

    //----- Half Length Postprocessing -----------------------------------------

    tr = x[0].re;
    ti = x[0].im;

    x[0].re = tr + ti;
    x[0].im = 0;

    x[nfft].re = tr - ti;
    x[nfft].im = 0;

    x[nfft/2].im = -x[nfft/2].im;

    for (i=1; i<nfft/2; i++)
    {
        j = nfft - i;

        rs = (x[i].re + x[j].re) * 0.5f;
        rd = (x[j].re - x[i].re) * 0.5f;
        is = (x[i].im + x[j].im) * 0.5f;
        id = (x[i].im - x[j].im) * 0.5f;

        ci = (float) (cos2table[i << nstride]);
        cj = (float) (cos2table[(nfft/2-i) << nstride]);

        rp = is * ci + rd * cj;
        ip = rd * ci - is * cj;

        x[i].re = (rp + rs);
        x[j].re = (rs - rp);

        x[i].im = (ip + id);
        x[j].im = (ip - id);
    }
}

//------------------------------------------------------------------------------

void irfft(complex_float32 *spectrum, float *output, int n)
{
    int i, j, k, nfft, nstride;
    float t0, tn, rs, is, rd, id, rp, ip, ci, cj, norm;
    complex_float32 *x;
    complex_float64 *w;
    double *cos2table;

    nfft = n/2;

    if (nfft == 0 || spectrum == NULL || output == NULL)
        return;

    if (ilog2(n) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", n);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(n);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(n);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2*k, 2*nfft);
        return;
    }

    x = (complex_float32 *) output;
    w = table.twiddle_factor;
    cos2table = table.cos_half;

    //----- Half Length Preprocessing ------------------------------------------

    t0 = spectrum[0].re;
    tn = spectrum[nfft].re;

    x[0].re = (t0 + tn);
    x[0].im = (t0 - tn);

    x[nfft/2].re = spectrum[nfft/2].re * 2;
    x[nfft/2].im = -spectrum[nfft/2].im * 2;

    for (i=1; i<nfft/2; i++)
    {
        j = nfft - i;

        rs = (spectrum[i].re + spectrum[j].re);
        rd = (spectrum[i].re - spectrum[j].re);
        is = (spectrum[i].im + spectrum[j].im);
        id = (spectrum[i].im - spectrum[j].im);

        ci = (float) (cos2table[i << nstride]);
        cj = (float) (cos2table[(nfft/2-i) << nstride]);

        rp = is * ci + rd * cj;
        ip = rd * ci - is * cj;

        x[i].re = rp + rs;
        x[j].re = rs - rp;

        x[i].im = ip - id;
        x[j].im = ip + id;
    }

    fft_core(x, w, nstride, nfft);

    norm = 1 / (float) n;

    for (i=0; i <n/2; i++)
    {
        output[2*i + 0] = x[i].re * norm;
        output[2*i + 1] = x[i].im * norm;
    }
}

//------------------------------------------------------------------------------

void cfft(complex_float32 *x, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    if (ilog2(nfft) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", nfft);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(nfft);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(nfft);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2 * k, 2 * nfft);
        return;
    }

    fft_core(x, table.twiddle_factor, nstride, nfft);
}

//------------------------------------------------------------------------------

void icfft(complex_float32 *x, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    if (ilog2(nfft) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", nfft);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(nfft);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(nfft);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2 * k, 2 * nfft);
        return;
    }

    for (i=0; i<nfft; i++)
        x[i].im = -x[i].im;

    fft_core(x, table.twiddle_factor, nstride, nfft);

    for (i=0; i<nfft; i++)
        x[i].im = -x[i].im;
}

//------------------------------------------------------------------------------

void magnitude(complex_float32 *input, float *result, int n)
{
    int i;

    for (i=0; i<n; i++)
    {
        result[i] = (input[i].re * input[i].re) + (input[i].im * input[i].im);
    }
}

//------------------------------------------------------------------------------

void magnitude_db(complex_float32 *input, float *result, int n)
{
    int i;
    double mag;

    for (i=0; i<n; i++)
    {
        mag = (input[i].re * input[i].re) + (input[i].im * input[i].im);

        if (mag < 1e-40)        // lower limit -400 dB
            mag = 1e-40;

        result[i] = (float) (10.*log10(mag));
    }
}

//------------------------------------------------------------------------------

void phase_rad(complex_float32 *input, float *result, int n)
{
    int i;

    for (i=0; i<n; i++)
    {
        result[i] = (float) atan2(input[i].im, input[i].re);
    }
}

//------------------------------------------------------------------------------

void magnitude_double(complex_float64 *input, double *result, int n)
{
    int i;

    for (i=0; i<n; i++)
    {
        result[i] = (input[i].re * input[i].re) + (input[i].im * input[i].im);
    }
}

//------------------------------------------------------------------------------

void magnitude_db_double(complex_float64 *input, double *result, int n)
{
    int i;
    double mag;

    for (i=0; i<n; i++)
    {
        mag = (input[i].re * input[i].re) + (input[i].im * input[i].im);

        if (mag < 1e-40)        // lower limit -400 dB
            mag = 1e-40;

        result[i] = 10.*log10(mag);
    }
}
//------------------------------------------------------------------------------

void phase_rad_double(complex_float64 *input, double *result, int n)
{
    int i;

    for (i=0; i<n; i++)
    {
        result[i]= atan2(input[i].im, input[i].re);
    }
}

//------------------------------------------------------------------------------

void fft_core_double(complex_float64 *x, complex_float64 *w, int nstride, int nfft)
{
    int i, ig, j, k, l, ngroups, nbutterflies;

    complex_float64 ctemp;

    //----- Bit Reverse Section ------------------------------------------------

    j = 0;
    for (i=0; i<nfft-1; i++)
    {
        if (i<j)
        {
            ctemp.re = x[j].re;
            ctemp.im = x[j].im;

            x[j].re = x[i].re;
            x[j].im = x[i].im;

            x[i].re = ctemp.re;
            x[i].im = ctemp.im;
        }

        k = nfft / 2;
        while (k <= j)
        {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    //----- First Stage --------------------------------------------------------

    i = 0;
    j = 1;
    ngroups = nfft >> 1;

    for (ig=0; ig<ngroups; ig++)
    {
        ctemp.re = x[j].re;
        ctemp.im = x[j].im;

        x[j].re = x[i].re - ctemp.re;
        x[j].im = x[i].im - ctemp.im;

        x[i].re = x[i].re + ctemp.re;
        x[i].im = x[i].im + ctemp.im;

        i = i + 2;
        j = j + 2;
    }

    //----- Final log2(nfft)-1 Stages ------------------------------------------

    nbutterflies = 2;

    for (ngroups = nfft >> 2; ngroups > 0; ngroups >>= 1)
    {
        i = 0;
        j = nbutterflies;

        for (ig=0; ig<ngroups; ig++)
        {
            k = 0;
            for (l=0; l<nbutterflies; l++)
            {
                // Butterfly
                ctemp.re = x[j].re * w[k].re - x[j].im * w[k].im;
                ctemp.im = x[j].im * w[k].re + x[j].re * w[k].im;

                x[j].re = x[i].re - ctemp.re;
                x[j].im = x[i].im - ctemp.im;

                x[i].re = x[i].re + ctemp.re;
                x[i].im = x[i].im + ctemp.im;

                k += (ngroups << nstride);
                i++;
                j++;
            }

            i = i + nbutterflies;
            j = j + nbutterflies;
        }

        nbutterflies <<= 1;
    }
}

//------------------------------------------------------------------------------

void rfft_double(double *input, complex_float64 *spectrum, int n)
{
    int i, j, k, nfft, nstride;
    double tr, ti, rs, is, rd, id, rp, ip, ci, cj;
    complex_float64 *x, *w;
    double *cos2table;

    nfft = n/2;

    if (nfft == 0 || input == NULL || spectrum == NULL)
        return;

    if (ilog2(n) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", n);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(n);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(n);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2*k, 2*nfft);
        return;
    }

    // copy memory if not in-place
    if (input != (double *)spectrum)
    {
        for (i=0; i <n/2; i++)
        {
            spectrum[i].re = input[2*i + 0];
            spectrum[i].im = input[2*i + 1];
        }
    }

    x = spectrum;
    w = table.twiddle_factor;
    cos2table = table.cos_half;

    fft_core_double(x, w, nstride, nfft);

    //----- Half Length Postprocessing -----------------------------------------

    tr = x[0].re;
    ti = x[0].im;

    x[0].re = tr + ti;
    x[0].im = 0;

    x[nfft].re = tr - ti;
    x[nfft].im = 0;

    x[nfft/2].im = -x[nfft/2].im;

    for (i=1; i<nfft/2; i++)
    {
        j = nfft - i;

        rs = (x[i].re + x[j].re) * 0.5;
        rd = (x[j].re - x[i].re) * 0.5;
        is = (x[i].im + x[j].im) * 0.5;
        id = (x[i].im - x[j].im) * 0.5;

        ci = cos2table[i << nstride];
        cj = cos2table[(nfft/2-i) << nstride];

        rp = is * ci + rd * cj;
        ip = rd * ci - is * cj;

        x[i].re = (rp + rs);
        x[j].re = (rs - rp);

        x[i].im = (ip + id);
        x[j].im = (ip - id);
    }
}

//------------------------------------------------------------------------------

void irfft_double(complex_float64 *spectrum, double *output, int n)
{
    int i, j, k, nfft, nstride;
    double t0, tn, rs, is, rd, id, rp, ip, ci, cj, norm;
    complex_float64 *x, *w;
    double *cos2table;

    nfft = n/2;

    if (nfft == 0 || spectrum == NULL || output == NULL)
        return;

    if (ilog2(n) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", n);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(n);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(n);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2*k, 2*nfft);
        return;
    }

    x = (complex_float64 *) output;
    w = table.twiddle_factor;
    cos2table = table.cos_half;

    //----- Half Length Preprocessing ------------------------------------------

    t0 = spectrum[0].re;
    tn = spectrum[nfft].re;

    x[0].re = (t0 + tn);
    x[0].im = (t0 - tn);

    x[nfft/2].re = spectrum[nfft/2].re * 2;
    x[nfft/2].im = -spectrum[nfft/2].im * 2;

    for (i=1; i<nfft/2; i++)
    {
        j = nfft - i;

        rs = (spectrum[i].re + spectrum[j].re);
        rd = (spectrum[i].re - spectrum[j].re);
        is = (spectrum[i].im + spectrum[j].im);
        id = (spectrum[i].im - spectrum[j].im);

        ci = cos2table[i << nstride];
        cj = cos2table[(nfft/2-i) << nstride];

        rp = is * ci + rd * cj;
        ip = rd * ci - is * cj;

        x[i].re = rp + rs;
        x[j].re = rs - rp;

        x[i].im = ip - id;
        x[j].im = ip + id;
    }

    fft_core_double(x, w, nstride, nfft);

    norm = 1 / (double) n;

    for (i=0; i<n/2; i++)
    {
        output[2 * i + 0] = x[i].re * norm;
        output[2 * i + 1] = x[i].im * norm;
    }
}

//------------------------------------------------------------------------------

void cfft_double(complex_float64 *x, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    if (ilog2(nfft) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", nfft);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(nfft);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(nfft);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2 * k, 2 * nfft);
        return;
    }

    fft_core_double(x, table.twiddle_factor, nstride, nfft);
}

//------------------------------------------------------------------------------

void icfft_double(complex_float64 *x, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    if (ilog2(nfft) == 0)
    {
        printf("Error: number of FFT bins must be a power of two (%d)\n", nfft);
        return;
    }

    // create new table if missing
    if (table.nfft == 0)
        set_twiddle_table(nfft);

    // table stride
    k = table.nfft;

    if (k < nfft)
    {
        // create new table if too small
        set_twiddle_table(nfft);
        nstride = 1;
    }
    else
    {
        // compute table stride
        nstride = -1;
        i = 0;
        while (k)
        {
            if (k == nfft)
            {
                nstride = i;
                break;
            }
            k = k >> 1;
            i++;
        }
    }

    if (nstride < 0)
    {
        printf("Error: invalid table size: %d (nfft: %d)\n", 2 * k, 2 * nfft);
        return;
    }

    for (i=0; i<nfft; i++)
        x[i].im = -x[i].im;

    fft_core_double(x, table.twiddle_factor, nstride, nfft);

    for (i=0; i<nfft; i++)
        x[i].im = -x[i].im;
}

//------------------------------------------------------------------------------

int ilog2(int iarg)
{
    int i, n;

    // integer base 2 logarithm
    n = 0;
    for (i=1; i<iarg; i*=2)
        n++;

    if (i != iarg)
        return 0;

    return n;
}

/*------------------------------License----------------------------------------*\
| Copyright (c) 1988-2012 Uwe Simmer                        					|
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
| DEALINGS IN THE SOFTWARE.  													|
|																				|
| https://opensource.org/licenses/mit-license.php								|
\*-----------------------------------------------------------------------------*/
