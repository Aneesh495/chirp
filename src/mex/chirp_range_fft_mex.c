/*
 * CHIRP RADAR TOOLBOX
 * MEX Kernel: Fast Range-FFT with Windowing
 *
 * Computes windowed range fast Fourier transform across all chirps and virtual
 * channels in a 3D radar data cube using an optimized Radix-2 Cooley-Tukey FFT.
 *
 * Function Signature:
 *   rangeCube = chirp_range_fft_mex(adcData, windowVector, nFft)
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Check if N is a power of 2 */
static inline int is_power_of_two(mwSize n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

/* Reverse bits of 32-bit integer up to log2(N) bits */
static inline mwSize bit_reverse(mwSize x, int log2n) {
    mwSize n = 0;
    for (int i = 0; i < log2n; ++i) {
        n = (n << 1) | (x & 1);
        x >>= 1;
    }
    return n;
}

/* In-place Radix-2 Decimation-In-Time Cooley-Tukey FFT */
static void radix2_fft(double *real, double *imag, mwSize n, int log2n,
                       const double *twiddleCos, const double *twiddleSin) {
    /* 1. Bit-reversal permutation */
    for (mwSize i = 0; i < n; ++i) {
        mwSize j = bit_reverse(i, log2n);
        if (j > i) {
            double tr = real[i]; real[i] = real[j]; real[j] = tr;
            double ti = imag[i]; imag[i] = imag[j]; imag[j] = ti;
        }
    }

    /* 2. Cooley-Tukey butterfly stages */
    for (int s = 1; s <= log2n; ++s) {
        mwSize m = (mwSize)1 << s;
        mwSize m2 = m >> 1;
        mwSize step = n >> s;

        for (mwSize k = 0; k < n; k += m) {
            for (mwSize j = 0; j < m2; ++j) {
                mwSize twIdx = j * step;
                double wr = twiddleCos[twIdx];
                double wi = twiddleSin[twIdx];

                mwSize idxB = k + j + m2;
                mwSize idxA = k + j;

                double br = real[idxB];
                double bi = imag[idxB];

                double tr = wr * br - wi * bi;
                double ti = wr * bi + wi * br;

                real[idxB] = real[idxA] - tr;
                imag[idxB] = imag[idxA] - ti;
                real[idxA] = real[idxA] + tr;
                imag[idxA] = imag[idxA] + ti;
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("chirp:range_fft:invalidInputs",
            "Usage: rangeCube = chirp_range_fft_mex(adcData, windowVector, nFft)");
    }

    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("chirp:range_fft:invalidType", "adcData must be double.");
    }

    mwSize numDims = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dims = mxGetDimensions(prhs[0]);
    mwSize nAdc = dims[0];
    mwSize nChirps = (numDims >= 2) ? dims[1] : 1;
    mwSize nVirt = (numDims >= 3) ? dims[2] : 1;

    mwSize nFft = (mwSize)mxGetScalar(prhs[2]);
    const double *win = mxGetPr(prhs[1]);
    mwSize winLen = mxGetNumberOfElements(prhs[1]);
    if (winLen < nAdc) {
        mexErrMsgIdAndTxt("chirp:range_fft:windowMismatch", "Window vector length must match ADC sample count.");
    }

    const double *prReal = mxGetPr(prhs[0]);
    const double *prImag = mxGetPi(prhs[0]);

    /* Allocate output array: [nFft x nChirps x nVirt] complex */
    mwSize outDims[3] = { nFft, nChirps, nVirt };
    plhs[0] = mxCreateNumericArray(3, outDims, mxDOUBLE_CLASS, mxCOMPLEX);
    double *outReal = mxGetPr(plhs[0]);
    double *outImag = mxGetPi(plhs[0]);

    mwSize totalCols = nChirps * nVirt;

    if (is_power_of_two(nFft)) {
        /* Compute log2(nFft) */
        int log2n = 0;
        mwSize temp = nFft;
        while (temp > 1) { temp >>= 1; log2n++; }

        /* Precompute twiddle factors: W_N^k = e^(-j * 2*pi*k / N) */
        double *twCos = (double *)malloc((nFft / 2) * sizeof(double));
        double *twSin = (double *)malloc((nFft / 2) * sizeof(double));
        if (!twCos || !twSin) {
            mexErrMsgIdAndTxt("chirp:range_fft:outOfMemory", "Twiddle factor allocation failed.");
        }

        double angleStep = -2.0 * M_PI / (double)nFft;
        for (mwSize i = 0; i < nFft / 2; ++i) {
            twCos[i] = cos(angleStep * (double)i);
            twSin[i] = sin(angleStep * (double)i);
        }

        /* Work buffers for each column */
        double *workR = (double *)malloc(nFft * sizeof(double));
        double *workI = (double *)malloc(nFft * sizeof(double));
        if (!workR || !workI) {
            free(twCos); free(twSin);
            mexErrMsgIdAndTxt("chirp:range_fft:outOfMemory", "Work buffer allocation failed.");
        }

        for (mwSize col = 0; col < totalCols; ++col) {
            const double *inR = prReal + col * nAdc;
            const double *inI = prImag ? (prImag + col * nAdc) : NULL;

            /* Copy windowed samples into work buffers and zero-pad to nFft */
            for (mwSize n = 0; n < nAdc; ++n) {
                double w = win[n];
                workR[n] = inR[n] * w;
                workI[n] = inI ? (inI[n] * w) : 0.0;
            }
            if (nFft > nAdc) {
                memset(workR + nAdc, 0, (nFft - nAdc) * sizeof(double));
                memset(workI + nAdc, 0, (nFft - nAdc) * sizeof(double));
            }

            /* In-place radix-2 FFT */
            radix2_fft(workR, workI, nFft, log2n, twCos, twSin);

            /* Copy to output destination */
            memcpy(outReal + col * nFft, workR, nFft * sizeof(double));
            memcpy(outImag + col * nFft, workI, nFft * sizeof(double));
        }

        free(twCos);
        free(twSin);
        free(workR);
        free(workI);
    } else {
        /* General DFT fallback for non-power-of-2 lengths */
        double *cosTable = (double *)malloc(nFft * nAdc * sizeof(double));
        double *sinTable = (double *)malloc(nFft * nAdc * sizeof(double));
        if (!cosTable || !sinTable) {
            mexErrMsgIdAndTxt("chirp:range_fft:outOfMemory", "Twiddle table allocation failed.");
        }

        double twoPiOverN = 2.0 * M_PI / (double)nFft;
        for (mwSize k = 0; k < nFft; ++k) {
            for (mwSize n = 0; n < nAdc; ++n) {
                double angle = twoPiOverN * (double)(k * n);
                cosTable[k + n * nFft] = cos(angle);
                sinTable[k + n * nFft] = -sin(angle);
            }
        }

        for (mwSize col = 0; col < totalCols; ++col) {
            const double *inR = prReal + col * nAdc;
            const double *inI = prImag ? (prImag + col * nAdc) : NULL;
            double *dstR = outReal + col * nFft;
            double *dstI = outImag + col * nFft;

            for (mwSize k = 0; k < nFft; ++k) {
                double sumR = 0.0;
                double sumI = 0.0;
                const double *cosPtr = cosTable + k;
                const double *sinPtr = sinTable + k;

                for (mwSize n = 0; n < nAdc; ++n) {
                    double w = win[n];
                    double xr = inR[n] * w;
                    double xi = inI ? (inI[n] * w) : 0.0;
                    double c = *cosPtr;
                    double s = *sinPtr;
                    cosPtr += nFft;
                    sinPtr += nFft;

                    sumR += (xr * c - xi * s);
                    sumI += (xr * s + xi * c);
                }
                dstR[k] = sumR;
                dstI[k] = sumI;
            }
        }

        free(cosTable);
        free(sinTable);
    }
}
