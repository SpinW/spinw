/*===========================================================
 * eig_omp.cpp - Diagonalises a stack of matrices using
 *               Lapack calls and OpenMP.
 *
 * [V,E]=eig_omp(H,'orth','sort');
 *     where H is n x n x l, 
 *           V is n x n x l 
 *       and E is n x l
 *
 * 'orth' and 'sort' are options. See eig_omp.m for full info
 *
 * This is a MEX-file for MATLAB.
 *
 * In addition the the mexFunction, the following 
 * subsidiary functions are defined:
 *
 * void fliplr()    - Flips matrix by columns
 * void quicksort() - Quicksort algorithm, returns indices
 * void sort()      - Sort and permutes eigenvalue/vectors
 * void orth()      - Orthogonalises degenerate eigenvecs
 * 
 * Original Author: M. D. Le  [duc.le@stfc.ac.uk]
 * $Revision$ ($Date$)
 *=========================================================*/

#include <cfloat>
#include <cstring>
#include <cmath>
#include <limits>
#include "mex.h"
#include "matrix.h"
#include "lapack.h"

#ifndef _OPENMP
void omp_set_num_threads(int nThreads) {};
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num()  0
#else
#include <omp.h>
#endif

// Define LAPACK functions depending on type
template <typename T>
void igeev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, T *a, const ptrdiff_t *lda, T *w, T *vl, const ptrdiff_t *ldvl,
    T *vr, const ptrdiff_t *ldvr, T *work, const ptrdiff_t *lwork, T *rwork, ptrdiff_t *info) {
    mexErrMsgIdAndTxt("eig_omp:wrongtype","This function is only defined for single and double floats.");
}
template <> void igeev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, float *a, const ptrdiff_t *lda, float *w, float *vl, const ptrdiff_t *ldvl,
    float *vr, const ptrdiff_t *ldvr, float *work, const ptrdiff_t *lwork, float *rwork, ptrdiff_t *info) {
    return cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}
template <> void igeev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, double *a, const ptrdiff_t *lda, double *w, double *vl, const ptrdiff_t *ldvl,
    double *vr, const ptrdiff_t *ldvr, double *work, const ptrdiff_t *lwork, double *rwork, ptrdiff_t *info) {
    return zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
}

template <typename T>
void geev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, T *a, const ptrdiff_t *lda, T *wr, T *wi, T *vl,
    const ptrdiff_t *ldvl, T *vr, const ptrdiff_t *ldvr, T *work, const ptrdiff_t *lwork, ptrdiff_t *info) {
    mexErrMsgIdAndTxt("eig_omp:wrongtype","This function is only defined for single and double floats.");
}
template <> void geev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, float *a, const ptrdiff_t *lda, float *wr, float *wi, float *vl,
    const ptrdiff_t *ldvl, float *vr, const ptrdiff_t *ldvr, float *work, const ptrdiff_t *lwork, ptrdiff_t *info) {
    return sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
template <> void geev(const char *jobvl, const char *jobvr, const ptrdiff_t *n, double *a, const ptrdiff_t *lda, double *wr, double *wi, double *vl,
    const ptrdiff_t *ldvl, double *vr, const ptrdiff_t *ldvr, double *work, const ptrdiff_t *lwork, ptrdiff_t *info) {
    return dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}

template <typename T>
void ievr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, T *a, const ptrdiff_t *lda,
    const T *vl, const T *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const T *abstol, ptrdiff_t *m, T *w, T *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, T *work, const ptrdiff_t *lwork, T *rwork, const ptrdiff_t *lrwork,
    ptrdiff_t *iwork, const ptrdiff_t *liwork, ptrdiff_t *info) {
    mexErrMsgIdAndTxt("eig_omp:wrongtype","This function is only defined for single and double floats.");
}
template <> void ievr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, float *a, const ptrdiff_t *lda,
    const float *vl, const float *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const float *abstol, ptrdiff_t *m, float *w, float *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, float *work, const ptrdiff_t *lwork, float *rwork, const ptrdiff_t *lrwork,
    ptrdiff_t *iwork, const ptrdiff_t *liwork, ptrdiff_t *info) {
    return cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}
template <> void ievr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, double *a, const ptrdiff_t *lda,
    const double *vl, const double *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const double *abstol, ptrdiff_t *m, double *w, double *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, double *work, const ptrdiff_t *lwork, double *rwork, const ptrdiff_t *lrwork,
    ptrdiff_t *iwork, const ptrdiff_t *liwork, ptrdiff_t *info) {
    return zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info);
}

template <typename T>
void evr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, T *a, const ptrdiff_t *lda,
    const T *vl, const T *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const T *abstol, ptrdiff_t *m, T *w, T *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, T *work, const ptrdiff_t *lwork, ptrdiff_t *iwork,
    const ptrdiff_t *liwork, ptrdiff_t *info) {
    mexErrMsgIdAndTxt("eig_omp:wrongtype","This function is only defined for single and double floats.");
}
template <> void evr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, float *a, const ptrdiff_t *lda,
    const float *vl, const float *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const float *abstol, ptrdiff_t *m, float *w, float *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, float *work, const ptrdiff_t *lwork, ptrdiff_t *iwork,
    const ptrdiff_t *liwork, ptrdiff_t *info) {
    return ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
template <> void evr(const char *jobz, const char *range, const char *uplo, const ptrdiff_t *n, double *a, const ptrdiff_t *lda,
    const double *vl, const double *vu, const ptrdiff_t *il, const ptrdiff_t *iu, const double *abstol, ptrdiff_t *m, double *w, double *z,
    const ptrdiff_t *ldz, ptrdiff_t *isuppz, double *work, const ptrdiff_t *lwork, ptrdiff_t *iwork,
    const ptrdiff_t *liwork, ptrdiff_t *info) {
    return dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}

template <typename T>
void gesvd(const char *jobu, const char *jobvt, const ptrdiff_t *m, const ptrdiff_t *n, T *a, const ptrdiff_t *lda, T *s,
    T *u, const ptrdiff_t *ldu, T *vt, const ptrdiff_t *ldvt, T *work, const ptrdiff_t *lwork, T *rwork, ptrdiff_t *info) {
    mexErrMsgIdAndTxt("eig_omp:wrongtype","This function is only defined for single and double floats.");
}
template <> void gesvd(const char *jobu, const char *jobvt, const ptrdiff_t *m, const ptrdiff_t *n, float *a, const ptrdiff_t *lda, float *s,
    float *u, const ptrdiff_t *ldu, float *vt, const ptrdiff_t *ldvt, float *work, const ptrdiff_t *lwork, float *rwork, ptrdiff_t *info) {
    return cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}
template <> void gesvd(const char *jobu, const char *jobvt, const ptrdiff_t *m, const ptrdiff_t *n, double *a, const ptrdiff_t *lda, double *s,
    double *u, const ptrdiff_t *ldu, double *vt, const ptrdiff_t *ldvt, double *work, const ptrdiff_t *lwork, double *rwork, ptrdiff_t *info) {
    return zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
}

// Flips a (column-major) matrix by columns, like the matlab function.
template <typename T>
void fliplr(T *M, mwSignedIndex m, mwSignedIndex n, T *vec, bool isreal)
{
    int ii;
    T val, *p0r, *p1r, *p0i, *p1i;
    // Just a row vector, reverse order of values.
    if(m==1) {
        if(isreal) {
            p0r = M; p1r = M+n-1;
            for(ii=0; ii<(n/2); ii++) {
                val = *p1r;
                *p1r-- = *p0r;
                *p0r++ = val;
            }
        }
        else {
            p0r = M; p1r = M+n*2-2;
            p0i = M+1; p1i = M+n*2-1;
            for(ii=0; ii<(n/2); ii++) {
                val = *p1r; *p1r-- = *p0r; *p0r++ = val;
                val = *p1i; *p1i-- = *p0i; *p0i++ = val;
            }
        }
    }
    // Actual matrix - assume column major
    else {
        if(isreal) {
            size_t msz = m*sizeof(T);
            for(ii=0; ii<(n/2); ii++) {
                memcpy(vec, M+(n-ii-1)*n, msz);
                memcpy(M+(n-ii-1)*n, M+ii*n, msz);
                memcpy(M+ii*n, vec, msz);
            }
        }
        else {
            mwSignedIndex n2 = n*2;
            size_t m2sz = 2*m*sizeof(T);
            for(ii=0; ii<(n/2); ii++) {
                memcpy(vec, M+(n-ii-1)*n2, m2sz);
                memcpy(M+(n-ii-1)*n2, M+ii*n2, m2sz);
                memcpy(M+ii*n2, vec, m2sz);
            }
        }
    }
}

// Quicksort modified from public domain implementation by Darel Rex Finley.
//      http://alienryderflex.com/quicksort/
template <typename T>
void quicksort(int *id, T *val, mwSignedIndex elements)
{
    T piv;
    int i=0, L, R, C, swap;
    int beg[300], end[300];
    beg[0]=0; end[0]=(int)elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R)
        {
            piv=val[id[L]]; C=id[L];
            while (L<R)
            {
                while (val[id[R]]>=piv && L<R) R--; if (L<R) id[L++]=id[R];
                while (val[id[L]]<=piv && L<R) L++; if (L<R) id[R--]=id[L];
            }
            id[L]=C; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
            if (i > 300) {
                mexErrMsgIdAndTxt("eig_omp:qsortoverflow", "Quicksort ran out of memory");
            }
            if (end[i]-beg[i]>end[i-1]-beg[i-1])
            {
                swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
                swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
            }
        }
        else i--;
    }
}

// This is called for real general matrices - complex conjugate eigenvectors stored in consecutive columns
template <typename T>
void sort(mwSignedIndex m, T *Dr, T *Di, T *V, T *work, int sort_type)
{
    int *id, ii, jj;
    T *val, *pD, *pDi, *pV; 
    T abstol = sqrt(std::numeric_limits<T>::epsilon());
    size_t msz;

    // Assume(!) that workspace is max[3*m+601,m*(m+3)] large.
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Probably should modify quicksort to allow descending sort
    if(sort_type==-1) {
        val = work+m+1;
        // For now just feed quicksort the negative of the eigenvalues...
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = -Dr[ii];
    }
    else 
        val = Dr;

    // Does the sorting using QuickSort (probably not the quickest algorithm, but fast and simple...)
    quicksort(id, val, (int)m);

    // Now permute the eigenvalues and eigenvectors 
    //   this is horendously bloated... maybe try Fich et al. 1995 
    //   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.29.2256
    msz = m*sizeof(T);
    pD = work+m; memcpy(pD, Dr, msz);
    pDi = work+2*m; memcpy(pDi, Di, msz);
    if(V != NULL) {
        pV = work+3*m; memcpy(pV, V, m*m*sizeof(T)); }
    for(ii=0; ii<(int)m; ii++) {
        // Take care to preserve the order of the eigenvectors (real parts first)
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol && Di[id[ii]]<abstol) {
                jj = id[ii+1];
                id[ii+1] = id[ii];
                id[ii] = jj;
            }
        pD[ii] = Dr[id[ii]];
        pDi[ii] = Di[id[ii]];
        if(V != NULL) {
            memcpy(&pV[ii*m], &V[id[ii]*m], msz); }
    }
    memcpy(Dr, pD, msz);
    memcpy(Di, pDi, msz);
    if(V != NULL) {
        memcpy(V, pV, m*m*sizeof(T)); }
}

// This is called for complex general matrices - both eigenvalues and eigenvectors are actually complex
template <typename T>
void sort(mwSignedIndex m, T *D, T *V, T *work, int sort_type)
{
    int *id, ii;
    T *val, *pD, *pV;
    size_t msz;

    // Assume(!) that workspace is max[3*m+601,2*m*(m+3)] large.
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Populate the sort variable as the real part of the eigenvalues
    val = work+m+1;
    if(sort_type==-1) {
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = -D[ii*2];
    }
    else {
        for(ii=0; ii<(int)m; ii++) 
            val[ii] = D[ii*2];
    }

    // Does the sorting using QuickSort
    quicksort(id, val, (int)m);

    // Now permute the eigenvalues and eigenvectors - in future look at using Fich et al. 1995
    msz = 2*m*sizeof(T);
    pD = work+m; memcpy(pD, D, msz);
    if(V != NULL) {
        pV = work+3*m; memcpy(pV, V, 2*m*m*sizeof(T)); }
    for(ii=0; ii<(int)m; ii++) {
        pD[ii*2] = D[id[ii]*2];
        pD[ii*2+1] = D[id[ii]*2+1];
        if(V != NULL) {
            memcpy(&pV[ii*2*m], &V[id[ii]*2*m], msz); }
    }
    memcpy(D, pD, msz);
    if(V != NULL) {
        memcpy(V, pV, 2*m*m*sizeof(T)); }
}

// This is called for real general matrices - complex conjugate eigenvectors stored in consecutive columns
template <typename T>
int orth(mwSignedIndex m, T *Dr, T *Di, T *Vr, T *Vi, T *work, bool isreal)
{
    int *id, ii, jj, kk, nn;
    mwSignedIndex n, info;
    T *pVz, *pS;
    T abstol = sqrt(std::numeric_limits<T>::epsilon());
    char jobu = 'O';
    char jobvt = 'N';
    mwSignedIndex lwork = 3*m;
    T *zwork = work + (2*m*(m+7)) - 2*3*m;
    T *rwork = zwork - 5*m;

    // Assume workspace is size [m*(m+7)]*sizeof(complex)
    id = (int*)work;
    for(ii=0; ii<(int)m; ii++) 
        id[ii] = ii;

    // Use QuickSort to get indices of eigenvalues in ascending order by real part
    quicksort(id, Dr, m);

    // Check for warning
//  for(ii=0; ii<(int)m; ii++)
//      if(fabs(Di[ii])>abstol) {
//          mexWarnMsgIdAndTxt("eig_omp:notdefinite","matrix contains complex eigenvalues, othogonalisation will not be accurate");
//          break;
//      }

    pVz = work+2*m;
    for(ii=0; ii<(int)m; ii++) {
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol) {
                if(Di[id[ii]]<-abstol) {
                    jj = id[ii+1];
                    id[ii+1] = id[ii];
                    id[ii] = jj;
                }
                kk = ii;
                while(true)
                    if(fabs(Dr[id[ii]]-Dr[id[++kk]])>abstol || kk >= (int)(m - 1))
                        break;
                // Populates a complex matrix with the degenerate eigenvectors, taking care about interleaving
                n = 0;
                for(jj=ii; jj<kk; jj++) {
                    // Conjugate pairs appear consecutively with positive imaginary part first
                    if(Di[id[jj]]>abstol) {
                        for(nn=0; nn<m; nn++) {
                            pVz[n*2*m+nn*2] = Vr[id[jj]*m+nn];
                            pVz[n*2*m+nn*2+1] = Vr[id[jj+1]*m+nn];
                            pVz[(n+1)*2*m+nn*2] = Vr[id[jj]*m+nn];
                            pVz[(n+1)*2*m+nn*2+1] = -Vr[id[jj+1]*m+nn];
                        }
                        jj++;
                        n++;
                    }
                    else
                        for(nn=0; nn<m; nn++) 
                            pVz[n*2*m+nn*2] = Vr[id[jj]*m+nn];
                    n++;
                }
                // Does the singular value decomposition to get the orthogonal basis and singular values
                pS = work + 2*(n+1)*m;
                gesvd(&jobu, &jobvt, &m, &n, pVz, &m, pS, pVz, &m, pVz, &m, zwork, &lwork, rwork, &info);
                // Checks that number of singular values == n
                for(nn=0; nn<n; nn++)
                    if(fabs(pS[nn])<0)
                        break;
                if((n-nn)>1)
                    return 1;
                // Re-interleave to Matlab output complex matrix format
                n = 0;
                for(jj=ii; jj<kk; jj++) {
                    for(nn=0; nn<m; nn++) {
                        Vr[id[jj]*m+nn] = pVz[n*2*m+nn*2];
                        Vi[id[jj]*m+nn] = pVz[n*2*m+nn*2+1];
                    }
                    n++;
                }
                ii = jj-1;
            }
    }
    return 0;
}

// This is called for complex general matrices - both eigenvalues and eigenvectors are actually complex
template <typename T>
int orth(mwSignedIndex m, T *D, T *V, T *work, bool isreal)
{
    int *id, ii, jj, kk, nn;
    mwSignedIndex n, info;
    T *Dr, *pVz, *pS;
    T abstol = sqrt(std::numeric_limits<T>::epsilon());
    char jobu = 'O';
    char jobvt = 'N';
    mwSignedIndex lwork = 3*m;
    T *zwork = work + (2*m*(m+10)) - 2*3*m;
    T *rwork = zwork - 2*5*m;
    size_t msz = 2*m*sizeof(T);

    // Assume workspace is size [m*(m+10)]*sizeof(complex)
    id = (int*)work;
    Dr = work + m;
    for(ii=0; ii<(int)m; ii++) {
        id[ii] = ii;
        Dr[ii] = D[2*ii];
    }

    // Use QuickSort to get indices of eigenvalues in ascending order by real part
    quicksort(id, Dr, m);

    pVz = work+2*m;
    for(ii=0; ii<(int)m; ii++) {
        if((ii+1)<(int)m)
            if(fabs(Dr[id[ii]]-Dr[id[ii+1]])<abstol) {
                kk = ii;
                while(true)
                    if(fabs(Dr[id[ii]]-Dr[id[++kk]])>abstol || kk>=(int)(m-1)) 
                        break;
                // Populates a complex matrix with the degenerate eigenvectors, taking care about interleaving
                n = 0;
                for(jj=ii; jj<kk; jj++)
                    memcpy(&pVz[(n++)*2*m], &V[id[jj]*2*m], msz);
                // Does the singular value decomposition to get the orthogonal basis and singular values
                pS = work + 2*m + 2*(n+1)*m;
                gesvd(&jobu, &jobvt, &m, &n, pVz, &m, pS, pVz, &m, pVz, &m, zwork, &lwork, rwork, &info);
                // Checks that number of singular values == n
                for(nn=0; nn<n; nn++)
                    if(fabs(pS[nn])<0)
                        break;
                if((n-nn)>1)
                    return 1;
                // Permutes back
                n = 0;
                for(jj=ii; jj<kk; jj++)
                    memcpy(&V[id[jj]*2*m], &pVz[(n++)*2*m], msz);
                ii = jj;
            }
    }
    return 0;
}

template <typename T>
void check_sym(bool *issym, const mxArray *mat, T tolsymm, mwSignedIndex m, int nthread, int *blkid, bool is_complex)
{
#pragma omp parallel default(none) shared(issym, mat) firstprivate(nthread, m, tolsymm, blkid, is_complex)
    if(is_complex) {
        T *A = (T*)mxGetData(mat);
        T *Ai = (T*)mxGetImagData(mat);
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            for(int ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                for(int ii=0; ii<m; ii++) {
                    // Just need to iterate over upper triangle
                    for(int jj=ii; jj<m; jj++) {
                        if(fabs( *(A+ii+jj*m) - *(A+jj+ii*m) ) > tolsymm
                         ||fabs( *(Ai+ii+jj*m) + *(Ai+jj+ii*m) ) > tolsymm) {
                            issym[ib] = false;
                            break;
                        }
                    }
                    if(!issym[ib])
                        break;
                }
                A += m*m;
                Ai += m*m;
            }
        }
    }
    else {
        T *A = (T*)mxGetData(mat);
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            for(int ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                for(int ii=0; ii<m; ii++) {
                    for(int jj=ii; jj<m; jj++) {
                        if(fabs( *(A+ii+jj*m) - *(A+jj+ii*m) ) > tolsymm) {
                            issym[ib] = false;
                            break;
                        }
                    }
                    if(!issym[ib])
                        break;
                }
                A += m*m;
            }
        }
    }
}

template <typename T>
int do_loop(T *mat, T *mat_i, mxArray *plhs[], int nthread, mwSignedIndex m, int nlhs, mwSignedIndex nd,
    const int *blkid, char jobz, bool anynonsym, const bool *issym, bool do_orth, int do_sort, bool is_complex);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex m, n, nd;
    const size_t *dims;
    int ib, nblock, nb, do_sort=0;
    int *blkid, valint;
    char jobz, *parstr, *valstr;
    bool *issym, anynonsym=false, do_orth=false, do_Colpa=false;
    bool is_single;
    int nthread = omp_get_max_threads();
    int err_code = 0;
//  mexPrintf("Number of threads = %d\n",nthread);

    // Checks inputs
    if(nrhs<1) {
        mexErrMsgIdAndTxt("eig_omp:nargin","Number of input argument must be at least 1.");
    }
    if(mxIsDouble(prhs[0])) {
        is_single = false;
    } else if(mxIsSingle(prhs[0])) {
        is_single = true;
    } else {
        mexErrMsgIdAndTxt("eig_omp:notfloat","Input matrix must be a float array.");
    }
    nd = mxGetNumberOfDimensions(prhs[0]);
    if(nd<2 || nd>3) {
        mexErrMsgIdAndTxt("eig_omp:wrongdims","Only 2D or 3D arrays are supported.");
    }
    dims = mxGetDimensions(prhs[0]);
    m = dims[0];
    n = dims[1];
    if(m!=n) {
        mexErrMsgIdAndTxt("eig_omp:notsquare","Input matrix is not square.");
    }
    if(nd==3)
        nblock = (int)dims[2];
    else
        nblock = 1;
    // Checks for optional arguments
    for(ib=1; ib<nrhs; ib++)
    {
        if(mxIsChar(prhs[ib])) {
            parstr = mxArrayToString(prhs[ib]);
            if(strcmp(parstr,"orth")==0) {
                if((nrhs-1)>ib) {
                    if(!mxIsChar(prhs[ib+1]) && (mxGetM(prhs[ib+1])*mxGetN(prhs[ib+1]))>0) {
                        if(mxIsLogical(prhs[ib+1]))
                            valint = *(int*)mxGetData(prhs[ib+1]);
                        else 
                            valint = (int)(*(double*)mxGetData(prhs[ib+1]));
                        do_orth = (valint!=0) ? true : false;
                        ib++;
                    }
                    else
                        do_orth = true;
                }
                else {
                    do_orth = true;
                }
            }
            else if(strcmp(parstr,"sort")==0) {
                if((nrhs-1)>ib) {
                    if(mxIsChar(prhs[ib+1])) {
                        valstr = mxArrayToString(prhs[ib+1]);
                        if(strcmp(valstr,"orth")==0)
                            continue;
                        else if(strcmp(valstr,"ascend")==0)
                            do_sort = 1;
                        else if(strcmp(valstr,"descend")==0)
                            do_sort = -1;
                        else 
                            mexErrMsgIdAndTxt("eig_omp:badsortarg","Arguments to 'sort' keyword must be either 'ascend' or 'descend'.");
                    }
                    else {
                        if(mxIsLogical(prhs[ib+1]))
                            valint = *(int*)mxGetData(prhs[ib+1]);
                        else 
                            valint = (int)(*(double*)mxGetData(prhs[ib+1]));
                        do_sort = (valint==-1||valint==0) ? valint : 1;
                    }
                    ib++;
                }
                else
                    do_sort = 1;
            }
            else if(strcmp(parstr,"Colpa")==0) {
                if((nrhs-1)>ib) {
                    if(!mxIsChar(prhs[ib+1]) && (mxGetM(prhs[ib+1])*mxGetN(prhs[ib+1]))>0) {
                        if(mxIsLogical(prhs[ib+1]))
                            valint = *(int*)mxGetData(prhs[ib+1]);
                        else 
                            valint = (int)(*(double*)mxGetData(prhs[ib+1]));
                        do_Colpa = (valint!=0) ? true : false;
                        ib++;
                    }
                    else
                        do_Colpa = true;
                }
                else {
                    do_Colpa = true;
                }
            }
        }
    }
    //mexPrintf("do_orth = %d; do_sort = %d\n",do_orth,do_sort);

    // More efficient to group blocks together to run in a single thread than to spawn one thread per matrix.
    if(nblock<nthread) {
        nthread = nblock;
    }
    blkid = new int[nthread+1];
    blkid[0] = 0;
    blkid[nthread] = nblock;
    nb = nblock/nthread;
    for(int nt=1; nt<nthread; nt++) {
       blkid[nt] = blkid[nt-1]+nb;
    }

    // Checks if all matrices are symmetric / hermitian.
    issym = new bool[nblock];
    // Initially assume all matrices symmetric 
    memset(issym,true,nblock*sizeof(bool));

    mxClassID classid;
    bool is_complex = mxIsComplex(prhs[0]);
    if(is_single) {
        // Tolerance on whether matrix is symmetric/hermitian
        check_sym(issym, prhs[0], (float)(FLT_EPSILON * 10.0), m, nthread, blkid, is_complex);
        classid = mxSINGLE_CLASS;
    } else {
        check_sym(issym, prhs[0], (double)(DBL_EPSILON * 10.0), m, nthread, blkid, is_complex);
        classid = mxDOUBLE_CLASS;
    }

    for(int ib=0; ib<nblock; ib++)
        if(!issym[ib]) {
            anynonsym = true;
            break;
        }

    // Some matrices are not symmetric, use [ZD]GEEV algorithm, get complex eigenvalues
    mxComplexity complexflag = anynonsym ? mxCOMPLEX : mxREAL;

    // Creates outputs
    if(nlhs<=1) {
        jobz = 'N';
        if(nd==2)
            plhs[0] = mxCreateNumericMatrix(m, 1, classid, complexflag);
        else
            plhs[0] = mxCreateNumericMatrix(m, nblock, classid, complexflag);
    }
    else {
        jobz = 'V';
        if(nd==2)
            plhs[1] = mxCreateNumericMatrix(m, m, classid, complexflag);
        else
            plhs[1] = mxCreateNumericMatrix(m, nblock, classid, complexflag);

        // If some matrices are not symmetric, will get complex conjugate eigenpairs
        complexflag = (mxIsComplex(prhs[0]) || anynonsym) ? mxCOMPLEX : mxREAL;
        if(nd==2)
            plhs[0] = mxCreateNumericMatrix(m, m, classid, complexflag);
        else
            plhs[0] = mxCreateNumericArray(3, dims, classid, complexflag);
    }
    //mexPrintf("IsComplex=%d, anynonsym=%d, nlhs=%d\n", mxIsComplex(prhs[0]), anynonsym, nlhs);
    //mexEvalString("drawnow;");

    void *i_part = is_complex ? mxGetImagData(prhs[0]) : NULL;
    if(is_single) {
        err_code = do_loop((float *)mxGetData(prhs[0]), (float *)i_part, plhs, nthread, m, nlhs, nd,
                           blkid, jobz, anynonsym, issym, do_orth, do_sort, is_complex);
    } else {
        err_code = do_loop((double *)mxGetData(prhs[0]), (double *)i_part, plhs, nthread, m, nlhs, nd,
                           blkid, jobz, anynonsym, issym, do_orth, do_sort, is_complex);
    }

    delete[]blkid; delete[]issym;
    if(err_code > 0)
        mexErrMsgIdAndTxt("eig_omp:defectivematrix","Eigenvectors of defective eigenvalues cannot be orthogonalised.");
}

template <typename T>
int do_loop(T *mat, T *mat_i, mxArray *plhs[], int nthread, mwSignedIndex m, int nlhs, mwSignedIndex nd,
    const int *blkid, char jobz, bool anynonsym, const bool *issym, bool do_orth, int do_sort, bool is_complex)
{
    int err_code = 0;
    T* lhs0 = (T*)mxGetData(plhs[0]);
    T* lhs1 = (T*)mxGetData(plhs[1]);
    T* ilhs0 = (T*)mxGetImagData(plhs[0]);
    T* ilhs1 = (T*)mxGetImagData(plhs[1]);
#pragma omp parallel default(none) shared(mat, mat_i, err_code, blkid, issym) \
    firstprivate(nthread, m, nlhs, nd, jobz, anynonsym, do_orth, do_sort, is_complex, lhs0, lhs1, ilhs0, ilhs1)
    {
#pragma omp for
        for(int nt=0; nt<nthread; nt++) {
            // These variables must be declared within the loop to make them local (and private)
            //   or we get memory errors.
            T *M, *E, *V=NULL, *D, *Di, *ptr_M, *ptr_Mi, *ptr_V, *ptr_Vi;
            mwSignedIndex m2, m22;
            size_t msz;
            char uplo = 'U';
            char range = 'A';
            char jobzn = 'N';
            T vu = std::numeric_limits<T>::max();
            T vl = -vu;
            mwSignedIndex il = 0, iu;
            T abstol = sqrt(std::numeric_limits<T>::epsilon());
            mwSignedIndex lda, ldz, numfind;
            mwSignedIndex info, lwork, liwork, lzwork;
            mwSignedIndex *isuppz, *iwork;
            T *work, *zwork;
            int ii, jj;
            lda = m;
            ldz = m;
            iu = m;
            m2 = m*m;
            m22 = 2*m2;
            msz = m2*sizeof(T);
            lwork = do_orth ? ( (26*m>(2*m*(m+10))) ? 26*m : (2*m*(m+10)) )  
                            : ( (26*m>(2*m*(m+3))) ? 26*m : (2*m*(m+3)) );
            liwork = 10*m;
            isuppz = new mwSignedIndex[2*m];
            work = new T[lwork];
            iwork = new mwSignedIndex[liwork];
            if(is_complex) {
                lzwork = 4*m;
                M = new T[m22];
                zwork = new T[lzwork*2];
                if(nlhs>1)
                    V = new T[m22];
            }
            else
                M = new T[m2];
            // The output of the _evr Lapack routines gives eigenvalues as vectors
            // If we want it as a diagonal matrix, need to use a temporary array...
            if(is_complex && anynonsym) {
                D = new T[2*m];
            }
            else if(nlhs>1 && nd==2) {
                D = new T[m];
                if(anynonsym)
                    Di = new T[m];
            }
            // Actual loop over individual matrices start here
            for(int ib=blkid[nt]; ib<blkid[nt+1]; ib++) {
                if(!(is_complex && anynonsym)) {
                    if(nlhs<=1)
                        D = lhs0 + ib*m;
                    else if(nd==3)
                        D = lhs1 + ib*m;
                }
                if(is_complex) {
                    ptr_M = mat + ib*m2;
                    ptr_Mi = mat_i + ib*m2;
                    // Interleaves complex matrices - Matlab stores complex matrix as an array of real
                    //   values followed by an array of imaginary values; Fortran (and C++ std::complex)
                    //   and hence Lapack stores it as arrays of pairs of values (real,imaginary).
                    for(ii=0; ii<m; ii++) {
                        for(jj=0; jj<m; jj++) {
                            M[ii*2+jj*2*m] = *ptr_M++;
                            M[ii*2+jj*2*m+1] = *ptr_Mi++;
                        }
                    }
                    // This does the actual computation - call using Matlab-included Intel MKL.
                    if(anynonsym) {
                        igeev(&jobz, &jobzn, &m, M, &lda, D, V, &ldz, V, &ldz, zwork, &lzwork, work, &info);
                        if(do_sort)
                            if(nlhs>1)
                                sort(m, D, V, work, do_sort);
                            else 
                                sort(m, D, (T*)NULL, work, do_sort);
                        if(nlhs>1 && do_orth)
                            if(orth(m, D, V, work, 0)==1) {
                                #pragma omp atomic
                                err_code++;
                                break;
                            }
                    }
                    else {
                        ievr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, zwork, &lzwork, work, &lwork, iwork, &liwork, &info);
                        // ZHEEVR outputs eigenvectors in ascending order by default.
                        if(do_sort==-1) {
                            fliplr(D,1,m,work,1);
                            if(nlhs>1) 
                                fliplr(V,m,m,work,0);
                        }
                    }
                    if(nlhs>1) {
                        ptr_V = lhs0 + ib*m2;
                        ptr_Vi = ilhs0 + ib*m2;
                        for(ii=0; ii<m; ii++) {
                            for(jj=0; jj<m; jj++) {
                                *ptr_V++ = V[ii*2*m+jj*2];
                                *ptr_Vi++ = -V[ii*2*m+jj*2+1];
                            }
                        }
                    }
                }
                else {
                    ptr_M = mat + ib*m2;
                    V = lhs0 + ib*m2;
                    if(anynonsym) {
                        if(nlhs<=1)
                            Di = ilhs0 + ib*m;
                        else if(nd==3)
                            Di = ilhs1 + ib*m;
                    }
                    // We cannot pass the right-side argument directly to Lapack because the routine
                    //   overwrites the input matrix to be diagonalised - so we create a copy
                    //   Matlab uses column-major like Fortran/Lapack (unlike C++) so can just memcpy
                    memcpy(M,ptr_M,msz);
                    // This does the actual computation - call using Matlab-included Intel MKL.
                    if(anynonsym) {
                        geev(&jobzn, &jobz, &m, M, &lda, D, Di, V, &ldz, V, &ldz, work, &lwork, &info);
                        if(do_sort)
                            if(nlhs>1)
                                sort(m, D, Di, V, work, do_sort);
                            else
                                sort(m, D, Di, (T*)NULL, work, do_sort);
                        if(nlhs>1 && do_orth)
                            if(orth(m, D, Di, V, lhs0+ib*m2, work, 1)==1) {
                                #pragma omp atomic
                                err_code++;
                                break;
                            }
                    }
                    else {
                        evr(&jobz, &range, &uplo, &m, M, &lda, &vl, &vu, &il, &iu, &abstol, &numfind, 
                                D, V, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
                        // DSYEVR outputs eigenvectors in ascending order by default.
                        if(do_sort==-1) {
                            fliplr(D,1,m,work,1);
                            if(nlhs>1) 
                                fliplr(V,m,m,work,1);
                        }
                    }
                    // Need to account for complex conjugate eigenvalue/vector pairs (imaginary parts 
                    //   of eigenvectors stored in consecutive columns)
                    if(nlhs>1 && !issym[ib] && !do_orth) {
                        ptr_V = lhs0 + ib*m2;
                        ptr_Vi = ilhs0 + ib*m2;
                        // Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue
                        //   having the positive imaginary part first.
                        for(ii=0; ii<m; ii++) {
                            if(*(Di+ii)>0.) {
                                for(jj=0; jj<m; jj++) {
                                    *(ptr_Vi+(ii+1)*m+jj) = -*(ptr_V+(ii+1)*m+jj);
                                    *(ptr_Vi+ii*m+jj) = *(ptr_V+(ii+1)*m+jj);
                                    *(ptr_V+(ii+1)*m+jj) = *(ptr_V+ii*m+jj);
                                }
                            }
                        }
                    }
                }
                // If the output needs eigenvalues as diagonal matrix, creates it.
                if(nlhs>1 && nd==2) { 
                    if(is_complex && anynonsym) {
                        E = lhs1 + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+2*ii);
                        E = ilhs1 + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+2*ii+1);
                    }
                    else {
                        E = lhs1 + ib*m2;
                        for(ii=0; ii<m; ii++) 
                            *(E+ii+ii*m) = *(D+ii);
                        if(anynonsym) {
                            E = ilhs1 + ib*m2; 
                            for(ii=0; ii<m; ii++) 
                                *(E+ii+ii*m) = *(Di+ii);
                        }
                    }
                }
                else if(is_complex && anynonsym) {
                    if(nlhs<=1)
                        E = lhs0 + ib*m;
                    else if(nd==3)
                        E = lhs1 + ib*m;
                    for(ii=0; ii<m; ii++) 
                        *(E+ii) = *(D+2*ii);
                    if(nlhs<=1)
                        E = ilhs0 + ib*m;
                    else if(nd==3)
                        E = ilhs1 + ib*m;
                    for(ii=0; ii<m; ii++) 
                        *(E+ii) = *(D+2*ii+1);
                }
                if(err_code!=0)
                    break;     // One of the threads has a defective matrix error - break loop here.
            }
            // Free memory...
            delete[]work; delete[]iwork; delete[]isuppz; delete[]M;
            if(is_complex)  {
                delete[]zwork;
                if(nlhs>1)
                    delete[]V;
            }
            if(is_complex && anynonsym) {
                delete[]D;
            }
            else if(nlhs>1 && nd==2) {
                delete[]D;
                if(anynonsym)
                    delete[]Di;
            }
            #ifndef _OPENMP
                if(err_code!=0)
                    break;
            #endif
        }
    }
    return err_code;
}
