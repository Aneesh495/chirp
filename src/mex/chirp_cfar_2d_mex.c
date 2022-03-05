/*
 * chirp_cfar_2d_mex.c
 * Native C MEX acceleration kernel for 2D Cell-Averaging Constant False Alarm Rate detection.
 *
 * Usage from MATLAB / Octave:
 *   [mask, thresh, nf] = chirp_cfar_2d_mex(rdmPower, pfa, trR, gdR, trV, gdV);
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 6) {
        mexErrMsgIdAndTxt("chirp:cfar:invalidInputs", 
            "Usage: [mask, thresh, nf] = chirp_cfar_2d_mex(rdmPower, pfa, trR, gdR, trV, gdV)");
    }

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("chirp:cfar:invalidType", "rdmPower must be a real double 2D matrix.");
    }

    mwSize nRows = mxGetM(prhs[0]);
    mwSize nCols = mxGetN(prhs[0]);
    const double *pwr = mxGetPr(prhs[0]);

    double pfa = mxGetScalar(prhs[1]);
    int trR = (int)mxGetScalar(prhs[2]);
    int gdR = (int)mxGetScalar(prhs[3]);
    int trV = (int)mxGetScalar(prhs[4]);
    int gdV = (int)mxGetScalar(prhs[5]);

    /* Allocate outputs */
    plhs[0] = mxCreateLogicalMatrix(nRows, nCols);
    mxLogical *mask = mxGetLogicals(plhs[0]);

    double *thresh = NULL;
    double *nf = NULL;

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
        thresh = mxGetPr(plhs[1]);
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(nRows, nCols, mxREAL);
        nf = mxGetPr(plhs[2]);
    }

    /* Total training cells count */
    int outerR = 2 * (trR + gdR) + 1;
    int outerV = 2 * (trV + gdV) + 1;
    int guardR = 2 * gdR + 1;
    int guardV = 2 * gdV + 1;
    int nTrain = (outerR * outerV) - (guardR * guardV);

    /* CA-CFAR alpha scale factor: alpha = N * (Pfa^(-1/N) - 1) */
    double alpha = (double)nTrain * (pow(pfa, -1.0 / (double)nTrain) - 1.0);

    /* Construct 2D Integral Image (1-based padded with zero borders) */
    mwSize intRows = nRows + 1;
    mwSize intCols = nCols + 1;
    double *intImg = (double *)calloc(intRows * intCols, sizeof(double));
    if (!intImg) {
        mexErrMsgIdAndTxt("chirp:cfar:outOfMemory", "Memory allocation failed for integral table.");
    }

    for (mwSize c = 0; c < nCols; ++c) {
        for (mwSize r = 0; r < nRows; ++r) {
            double val = pwr[r + c * nRows];
            double left = intImg[(r + 1) + c * intRows];
            double top = intImg[r + (c + 1) * intRows];
            double diag = intImg[r + c * intRows];
            intImg[(r + 1) + (c + 1) * intRows] = val + left + top - diag;
        }
    }

    int halfR = trR + gdR;
    int halfV = trV + gdV;

    /* Sliding window evaluation */
    for (int c = 0; c < (int)nCols; ++c) {
        int v1_out = (c - halfV < 0) ? 0 : c - halfV;
        int v2_out = (c + halfV >= (int)nCols) ? (int)nCols - 1 : c + halfV;
        int v1_in  = (c - gdV < 0) ? 0 : c - gdV;
        int v2_in  = (c + gdV >= (int)nCols) ? (int)nCols - 1 : c + gdV;

        for (int r = 0; r < (int)nRows; ++r) {
            int r1_out = (r - halfR < 0) ? 0 : r - halfR;
            int r2_out = (r + halfR >= (int)nRows) ? (int)nRows - 1 : r + halfR;
            int r1_in  = (r - gdR < 0) ? 0 : r - gdR;
            int r2_in  = (r + gdR >= (int)nRows) ? (int)nRows - 1 : r + gdR;

            /* Outer rectangle sum */
            double sumOut = intImg[(r2_out + 1) + (v2_out + 1) * intRows]
                          - intImg[r1_out + (v2_out + 1) * intRows]
                          - intImg[(r2_out + 1) + v1_out * intRows]
                          + intImg[r1_out + v1_out * intRows];

            /* Inner guard rectangle sum */
            double sumIn  = intImg[(r2_in + 1) + (v2_in + 1) * intRows]
                          - intImg[r1_in + (v2_in + 1) * intRows]
                          - intImg[(r2_in + 1) + v1_in * intRows]
                          + intImg[r1_in + v1_in * intRows];

            int areaOut = (r2_out - r1_out + 1) * (v2_out - v1_out + 1);
            int areaIn  = (r2_in - r1_in + 1) * (v2_in - v1_in + 1);
            int effectiveTrain = areaOut - areaIn;

            double localNf = (effectiveTrain > 0) ? ((sumOut - sumIn) / (double)effectiveTrain) : 1e-12;
            double localThresh = alpha * localNf;

            size_t idx = (size_t)r + (size_t)c * nRows;
            mask[idx] = (pwr[idx] > localThresh);

            if (thresh) thresh[idx] = localThresh;
            if (nf) nf[idx] = localNf;
        }
    }

    free(intImg);
}
