/*
 * chirp_backprojection_mex.c
 * Native C MEX acceleration kernel for time-domain SAR backprojection image formation.
 *
 * Usage from MATLAB / Octave:
 *   sarImage = chirp_backprojection_mex(rangeProfiles, rangeAxis, platformPos, gridX, gridY, groundZ, lambda);
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 7) {
        mexErrMsgIdAndTxt("chirp:bpa:invalidInputs", 
            "Usage: sarImage = chirp_backprojection_mex(rangeProfiles, rangeAxis, platformPos, gridX, gridY, groundZ, lambda)");
    }

    mwSize nRange = mxGetM(prhs[0]);
    mwSize nPulses = mxGetN(prhs[0]);
    const double *profR = mxGetPr(prhs[0]);
    const double *profI = mxGetPi(prhs[0]);

    const double *rAxis = mxGetPr(prhs[1]);
    double r0 = rAxis[0];
    double dr = (nRange > 1) ? (rAxis[1] - rAxis[0]) : 1.0;

    const double *platPos = mxGetPr(prhs[2]); /* 3 x nPulses */

    mwSize nX = mxGetNumberOfElements(prhs[3]);
    const double *gridX = mxGetPr(prhs[3]);

    mwSize nY = mxGetNumberOfElements(prhs[4]);
    const double *gridY = mxGetPr(prhs[4]);

    double groundZ = mxGetScalar(prhs[5]);
    double lambda = mxGetScalar(prhs[6]);
    double kWave = 4.0 * M_PI / lambda;

    /* Allocate output image: [nY x nX] complex double */
    plhs[0] = mxCreateDoubleMatrix(nY, nX, mxCOMPLEX);
    double *imgR = mxGetPr(plhs[0]);
    double *imgI = mxGetPi(plhs[0]);

    /* Precompute grid coordinates */
    for (mwSize u = 0; u < nPulses; ++u) {
        double px = platPos[0 + 3 * u];
        double py = platPos[1 + 3 * u];
        double pz = platPos[2 + 3 * u];

        const double *curProfR = profR + u * nRange;
        const double *curProfI = profI ? (profI + u * nRange) : NULL;

        for (mwSize ix = 0; ix < nX; ++ix) {
            double gx = gridX[ix];
            double dx = gx - px;

            for (mwSize iy = 0; iy < nY; ++iy) {
                double gy = gridY[iy];
                double dy = gy - py;
                double dz = groundZ - pz;

                double slantR = sqrt(dx * dx + dy * dy + dz * dz);
                double binPos = (slantR - r0) / dr;

                if (binPos >= 0.0 && binPos < (double)(nRange - 1)) {
                    int binFloor = (int)binPos;
                    double frac = binPos - (double)binFloor;

                    double s1_r = curProfR[binFloor];
                    double s1_i = curProfI ? curProfI[binFloor] : 0.0;
                    double s2_r = curProfR[binFloor + 1];
                    double s2_i = curProfI ? curProfI[binFloor + 1] : 0.0;

                    /* Linear interpolation */
                    double sig_r = (1.0 - frac) * s1_r + frac * s2_r;
                    double sig_i = (1.0 - frac) * s1_i + frac * s2_i;

                    /* Matched phase compensation: exp(j * kWave * slantR) */
                    double phase = kWave * slantR;
                    double cosP = cos(phase);
                    double sinP = sin(phase);

                    /* Complex multiplication: (sig_r + j*sig_i) * (cosP + j*sinP) */
                    double term_r = sig_r * cosP - sig_i * sinP;
                    double term_i = sig_r * sinP + sig_i * cosP;

                    size_t outIdx = iy + ix * nY;
                    imgR[outIdx] += term_r;
                    imgI[outIdx] += term_i;
                }
            }
        }
    }
}
