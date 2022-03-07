/*
 * chirp_range_fft_mex.c
 * Native C MEX acceleration kernel for fast-time windowing and Range-FFT.
 *
 * Usage from MATLAB / Octave:
 *   rangeProfiles = chirp_range_fft_mex(radarCube, windowWeights, nFft);
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("chirp:range_fft:invalidInputs", 
            "Usage: rangeProfiles = chirp_range_fft_mex(radarCube, windowWeights, nFft)");
    }

    const mwSize *dims = mxGetDimensions(prhs[0]);
    mwSize numDims = mxGetNumberOfDimensions(prhs[0]);
    
    mwSize nAdc = dims[0];
    mwSize nChirps = (numDims >= 2) ? dims[1] : 1;
    mwSize nVirt = (numDims >= 3) ? dims[2] : 1;

    mwSize nFft = (mwSize)mxGetScalar(prhs[2]);
    const double *win = mxGetPr(prhs[1]);

    const double *prReal = mxGetPr(prhs[0]);
    const double *prImag = mxGetPi(prhs[0]);

    /* Allocate output array: [nFft x nChirps x nVirt] complex */
    mwSize outDims[3] = { nFft, nChirps, nVirt };
    plhs[0] = mxCreateNumericArray(3, outDims, mxDOUBLE_CLASS, mxCOMPLEX);
    double *outReal = mxGetPr(plhs[0]);
    double *outImag = mxGetPi(plhs[0]);

    /* Precompute twiddle factors for DFT */
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
            sinTable[k + n * nFft] = -sin(angle); /* e^(-j * 2*pi*k*n / N) */
        }
    }

    /* Process all columns */
    mwSize totalCols = nChirps * nVirt;

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

                /* (xr + j*xi) * (c + j*s) = (xr*c - xi*s) + j*(xr*s + xi*c) */
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
