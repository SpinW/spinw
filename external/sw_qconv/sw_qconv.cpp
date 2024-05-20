#include "mex.h"
#include <stdexcept>
#include <complex>
#include <thread>
#include <vector>

template <typename T>
void loop(T *swOut, const mwSize *d1, const double stdG, const double *Qconv, const T *swConv, size_t i0, size_t i1) {
    for (size_t ii=i0; ii<i1; ii++) {
        double sumfG = 0.0;
        std::vector<double> fG(d1[1], 0.0);
        for (size_t jj=0; jj<d1[1]; jj++) {
            fG[jj] = exp(-pow((Qconv[jj] - Qconv[ii]) / stdG, 2) / 2);
            sumfG += fG[jj];
        }
        for (size_t jj=0; jj<d1[1]; jj++) {
            fG[jj] /= sumfG;
            for (size_t kk=0; kk<d1[0]; kk++) {
                swOut[kk+d1[0]*jj] += swConv[kk+d1[0]*ii] * fG[jj];
            }
        }
    }
}

template <typename T>
void do_calc(T *swOut, const mwSize *d1, const double stdG, const double *Qconv, const T *swConv, size_t nThread) {
    if (d1[1] > 10*nThread) {
        size_t nBlock = d1[1] / nThread;
        size_t i0 = 0, i1 = nBlock;
        std::vector<T*> swv(nThread);
        std::vector<std::thread> threads;
        for (size_t ii=0; ii<nThread; ii++) {
            swv[ii] = new T[d1[0] * d1[1]];
            threads.push_back(
                std::thread(loop<T>, std::ref(swv[ii]), std::ref(d1), stdG, std::ref(Qconv), std::ref(swConv), i0, i1)
            );
            i0 = i1;
            i1 += nBlock;
            if (i1 > d1[1] || ii == (nThread - 2)) {
                i1 = d1[1]; }
        }
        for (size_t ii=0; ii<nThread; ii++) {
            if (threads[ii].joinable()) {
                threads[ii].join(); }
            for (size_t jj=0; jj<d1[0]; jj++) {
                for (size_t kk=0; kk<d1[1]; kk++) {
                    swOut[jj+kk*d1[0]] += swv[ii][jj+kk*d1[0]];
                }
            }
            delete[](swv[ii]);
        }
    } else {
        loop<T>(swOut, d1, stdG, Qconv, swConv, 0, d1[1]);
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    if (nrhs < 3) {
        throw std::runtime_error("sw_qconv: Requires 3 arguments");
    }
    size_t nThreads = 8;
    if (nrhs == 4) {
        nThreads = (size_t)(*mxGetDoubles(prhs[3]));
    }
    if (mxIsComplex(prhs[1])) { throw std::runtime_error("Arg 2 is complex\n"); }
    if (mxIsComplex(prhs[2])) { throw std::runtime_error("Arg 3 is complex\n"); }
    const mwSize *d1 = mxGetDimensions(prhs[0]);
    const mwSize *d2 = mxGetDimensions(prhs[1]);
    if (d1[1] != d2[1]) { throw std::runtime_error("Arg 1 and 2 size mismatch\n"); }
    if (mxGetNumberOfElements(prhs[2]) > 1) { throw std::runtime_error("Arg 3 should be scalar\n"); }
    double *Qconv = mxGetDoubles(prhs[1]);
    double stdG = *(mxGetDoubles(prhs[2]));
    if (mxIsComplex(prhs[0])) {
        std::complex<double> *swConv = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));
        plhs[0] = mxCreateDoubleMatrix(d1[0], d1[1], mxCOMPLEX);
        std::complex<double> *swOut = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(plhs[0]));
        do_calc(swOut, d1, stdG, Qconv, swConv, nThreads);
    } else {
        double *swConv = mxGetDoubles(prhs[0]);
        plhs[0] = mxCreateDoubleMatrix(d1[0], d1[1], mxREAL);
        double *swOut = mxGetDoubles(plhs[0]);
        do_calc(swOut, d1, stdG, Qconv, swConv, nThreads);
    }
}
