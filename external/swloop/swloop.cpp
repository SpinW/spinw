#include <complex>
#include "mex.h"
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Eigenvalues"
#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <cmath>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
typedef std::complex<double> cd_t;
typedef std::tuple<std::complex<double>, Eigen::VectorXcd> ev_tuple_t;

struct pars {
    bool hermit;
    bool formfact;
    bool incomm;
    bool helical;
    bool bq;
    bool field;
    double omega_tol;
    size_t nMagExt;
    size_t nTwin;
    size_t nHkl;
    size_t nBond;
    size_t nBqBond;
    bool fastmode;
    bool neutron_output;
    size_t nformula;
    double nCell;
};

// In terms of the Toth & Lake paper ( https://arxiv.org/pdf/1402.6069.pdf ),
// ABCD = [A B; B' A'] in eq (25) and (26)
// ham_diag = [C 0; 0 C] in eq (25) and (26)
// zed = u^alpha in eq (9)

struct swinputs {
    const double *hklExt;                   // 3*nQ array of Q-points
    const std::complex<double> *ABCD;       // 3*nBond flattened array of non-q-dep part of H
    const double *idxAll;                   // indices into ABCD to compute H using accumarray
    const double *ham_diag;                 // the diagonal non-q-dep part of H
    const double *dR;                       // 3*nBond array of bond vectors
    const double *RR;                       // 3*nAtom array of magnetic ion coordinates
    const double *S0;                       // nAtom vector of magnetic ion spin lengths
    const std::complex<double> *zed;        // 3*nAtom array of moment direction vectors
    const double *FF;                       // nQ*nAtom array of magnetic form factor values
    const double *bqdR;                     // 3*nBond_bq array of biquadratic bond vectors
    const double *bqABCD;                   // 3*nBond_bq array of non-q-dep part biquadratic H
    const double *idxBq;                    // indices into bqABCD to compute H with accumarray
    const double *bq_ham_diag;              // diagonal part of the biquadratic Hamiltonian
    std::vector< const double *> ham_MF_v;  // The full Zeeman hamiltonian.
    const int *idx0;                        // Indices into hklExt for twined Q
    const double *n;                        // normal vector defining the rotating frame (IC calcs)
    const double *rotc;                     // 3*3*nTwin twin rotation matrices
    const double *hklA;                     // 3*nQ array of Q-points in Cartesian A^-1 units
};

struct swoutputs {
    std::complex<double> *omega;
    std::complex<double> *Sab;
    bool warn_posdef;
    bool warn_orth;
};

template <typename T>
using MatrixXT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
MatrixXT<T> accumarray(Eigen::MatrixXd ind,
                       Eigen::Array<T, Eigen::Dynamic, 1> data,
                       size_t sz1, size_t sz2=0) {
    // Note this version of accumarray is explicitly 2D
    if (sz2 == 0)
        sz2 = sz1;
    MatrixXT<T> res = MatrixXT<T>::Zero(sz1, sz2);
    for (int ii=0; ii<ind.rows(); ii++) {
        // Indices are direct from Matlab so use base-1 and are also doubles not ints
        res(static_cast<size_t>(ind(ii, 0)) - 1, static_cast<size_t>(ind(ii, 1)) - 1) += data(ii);
    }
    return res;
}

std::atomic_int err_flag(0);
std::vector<std::string> errmsgs = {
    std::string("all good"),
    std::string("swloop:notposdef: The input matrix is not positive definite."),
    std::string("swloop:cantdiag: Eigensolver could not converge")
};

void swcalc_fun(size_t i0, size_t i1, struct pars &params, struct swinputs &inputs, struct swoutputs &outputs) {
    // Setup
    size_t nHam = 2 * params.nMagExt;
    size_t nHam_sel = params.fastmode ? params.nMagExt : nHam;
    size_t len_abcd = 3 * params.nBond;
    // Assuming all data is column-major(!), so can be mapped directly to Eigen matrices
    // We're also assuming that the memory blocks are actually as big as defined here!
    // (E.g. the code which calls this needs to do some bounds checking)
    Eigen::Map<const Eigen::MatrixXd> hklExt(inputs.hklExt, 3, params.nHkl);
    Eigen::Map<const Eigen::VectorXcd> ABCD (inputs.ABCD, len_abcd);
    Eigen::Map<const Eigen::MatrixXd> idxAll (inputs.idxAll, len_abcd, 2);
    Eigen::Map<const Eigen::VectorXd> ham_diag (inputs.ham_diag, nHam);
    Eigen::Map<const RowMatrixXd> dR (inputs.dR, params.nBond, 3);
    Eigen::Map<const Eigen::MatrixXd> RR (inputs.RR, 3, params.nMagExt);
    Eigen::VectorXd sqrtS0 = Eigen::VectorXd(params.nMagExt);
    for (int ii=0; ii<sqrtS0.rows(); ii++) {
        sqrtS0(ii) =  sqrt(inputs.S0[ii] / 2.); }
    Eigen::Map<const Eigen::MatrixXcd> zed (inputs.zed, 3, params.nMagExt);
    Eigen::Map<const Eigen::MatrixXd> FF (inputs.FF, params.nMagExt, params.nHkl);
    Eigen::Map<const RowMatrixXd> bqdR (inputs.bqdR, params.nBqBond, 3);
    Eigen::Map<const Eigen::VectorXd> bqABCD (inputs.bqABCD, 3 * params.nBqBond);
    Eigen::Map<const Eigen::MatrixXd> idxBq (inputs.idxBq, 3 * params.nBqBond, 2);
    Eigen::Map<const Eigen::VectorXd> bq_ham_diag (inputs.bq_ham_diag, nHam);
    std::vector< Eigen::Map<const Eigen::MatrixXd> > ham_MF;
    if (params.field) {
        for (size_t ii=0; ii<params.nTwin; ii++) {
            ham_MF.push_back(Eigen::Map<const Eigen::MatrixXd>(inputs.ham_MF_v[ii], nHam, nHam));
        }
    }

    Eigen::MatrixXd gComm = Eigen::MatrixXd::Identity(nHam, nHam);
    for (size_t ii=params.nMagExt; ii<nHam; ii++) {
        gComm(ii, ii) = -1.; }
    Eigen::Map<const Eigen::VectorXi> idx0(inputs.idx0, params.nHkl);
    size_t nHklT = params.nHkl / params.nTwin;

    std::complex<double> *omega = outputs.omega;
    std::complex<double> *Sab_ptr, *Sperp_ptr;
    if (params.neutron_output) {
        Sab_ptr = new std::complex<double>[9];
        Sperp_ptr = outputs.Sab;
    } else {
        Sab_ptr = outputs.Sab;
    }

    size_t nHklI = params.nHkl / params.nTwin / 3;
    Eigen::Matrix3cd K1, K2, cK1, nx, K, qPerp;
    Eigen::Vector3d n;
    Eigen::Matrix3cd m1 = Eigen::Matrix3cd::Identity(3, 3);
    if (params.incomm) {
        // Defines the Rodrigues rotation matrix to rotate Sab back
        n << inputs.n[0], inputs.n[1], inputs.n[2];
        nx << 0, -n(2), n(1),
              n(2), 0, -n(0),
             -n(1), n(0), 0;
        K2 = n * n.adjoint().eval();
        K1 = 0.5 * (m1 - K2 - nx * std::complex<double>(0., 1.));
        cK1 = K1.conjugate().eval();
    }

    // This is the main loop
    for (size_t jj=i0; jj<i1; jj++) {
        if (jj % 10 == 0) {
            if (err_flag.load(std::memory_order_relaxed) > 0) {
                return; }
        }
        Eigen::MatrixXcd V(nHam, nHam_sel);
        Eigen::VectorXcd ExpF = exp((dR * hklExt(Eigen::all, jj)).array() * std::complex<double>(0., 1.));
        Eigen::MatrixXcd ham = accumarray<cd_t>(idxAll, ABCD.array() * ExpF.replicate(3,1).array(), nHam);
        ham += ham_diag.asDiagonal();
        if (params.bq) {
            Eigen::VectorXcd bqExp = exp((bqdR * hklExt(Eigen::all, jj)).array() * std::complex<double>(0., 1.));
            Eigen::MatrixXcd bqham = accumarray<cd_t>(idxBq, bqABCD.array() * bqExp.replicate(3,1).array(), nHam);
            ham += bqham;
            ham += bq_ham_diag.asDiagonal();
        }
        if (params.field) {
            ham += ham_MF[jj / nHklT];
        }
        ham = (ham + ham.adjoint().eval()) / 2.;
        if (params.hermit) {
            Eigen::LLT<Eigen::MatrixXcd> chol(ham);
            if (chol.info() == Eigen::NumericalIssue) {
                Eigen::VectorXcd eigvals = ham.eigenvalues();
                double tol0 = abs(eigvals.real().minCoeff()) * sqrt(nHam) * 4.;
                if (tol0 > params.omega_tol) {
                    err_flag.store(1, std::memory_order_relaxed);
                    return;
                }
                Eigen::MatrixXcd hamtol = ham + (Eigen::MatrixXcd::Identity(nHam, nHam) * tol0);
                chol.compute(hamtol);
                if (chol.info() == Eigen::NumericalIssue) {
                    chol.compute(ham + (Eigen::MatrixXcd::Identity(nHam, nHam) * params.omega_tol));
                    if (chol.info() == Eigen::NumericalIssue) {
                        err_flag.store(1, std::memory_order_relaxed);
                        return;
                    }
                }
                outputs.warn_posdef = true;
            }
            Eigen::MatrixXcd K = chol.matrixU();
            Eigen::MatrixXcd K2 = K * gComm * K.adjoint().eval();
            K2 = (K2 + K2.adjoint().eval()) / 2;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(K2);
            if (eig.info() != Eigen::Success) {
                err_flag.store(2, std::memory_order_relaxed);
                return;
            }
            std::vector<ev_tuple_t> evs;
            for (size_t ii=0; ii<nHam; ii++) {
                evs.push_back(ev_tuple_t(eig.eigenvalues()[ii], eig.eigenvectors().col(ii))); }
            std::sort(evs.begin(), evs.end(),
                [&](const ev_tuple_t &a, const ev_tuple_t &b) -> bool {
                    return std::real(std::get<0>(a)) > std::real(std::get<0>(b)); });
            Eigen::MatrixXcd U(nHam, nHam_sel);
            for (size_t ii=0; ii<nHam_sel; ii++) {
                // Note the following pointer arithmetic assumes column-major ordering
                // We assign the current eigenvalue to the element pointed to by omega (*omega=evs[ii])
                // then increment the pointer to point to the next element (omega++).
                // Since we do this q-point by q-point and the Matlab omega array has each
                // column having all the eigenvalues of a give Q-point, we need a column major layout.
                // If Matlab switches to a row-major layout, the unit tests will fail and we would
                // need to transpose the matrix after all these assignments.
                *(omega++) = std::get<0>(evs[ii]);
                U.col(ii) = std::get<1>(evs[ii]) * sqrt(gComm(ii,ii) * std::get<0>(evs[ii]));
            }
            V = K.inverse() * U;
        } else {
            Eigen::MatrixXcd gham = (gComm * ham).eval();
            // Add a small amount to the diagonal to ensure there are no degenerate levels,
            // so that the eigenvectors would be (quasi)orthogonal (instead of using eigorth()).
            double vd = 0.5 - (double)params.nMagExt;
            for (size_t ii=0; ii<nHam; ii++) { gham(ii, ii) += vd*1e-12; vd++; }
            Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eig(gham);
            if (eig.info() != Eigen::Success) {
                err_flag.store(2, std::memory_order_relaxed);
                return;
            }
            std::vector<ev_tuple_t> evs;
            for (size_t ii=0; ii<nHam; ii++) {
                evs.push_back(ev_tuple_t(eig.eigenvalues()[ii], eig.eigenvectors().col(ii))); }
            std::sort(evs.begin(), evs.end(),
                [&](const ev_tuple_t &a, const ev_tuple_t &b) -> bool {
                    return std::real(std::get<0>(a)) > std::real(std::get<0>(b)); });
            for (size_t ii=0; ii<nHam_sel; ii++) {
                *(omega++) = std::get<0>(evs[ii]);  // Note this pointer arithmetic assumes column-major ordering
                Eigen::ArrayXcd U = std::get<1>(evs[ii]);
                std::complex<double> vsum = (gComm.diagonal().array() * U.conjugate() * U).sum();
                V.col(ii) = U * sqrt(1.0 / vsum);
            }
        }
        Eigen::MatrixXcd zedExpF(nHam, 3);
        for (size_t j1=0; j1<params.nMagExt; j1++) {
            double ph = 0.;
            for (size_t i2=0; i2<3; i2++) {
                ph += hklExt(i2, idx0(jj)) * RR(i2, j1); }
            std::complex<double> expF = exp(std::complex<double>(0., -ph)) * sqrtS0(j1);
            if (params.formfact) {
                expF *= FF(j1, jj); }
            for (size_t i2=0; i2<3; i2++) {
                // Don't know why we have to use the reverse of the Matlab code here and in
                // V.tranpose() below instead of V.adjoint() - but otherwise get wrong intensities...
                zedExpF(j1, i2) = std::conj(zed(i2, j1)) * expF;
                zedExpF(j1 + params.nMagExt, i2) = zed(i2, j1) * expF;
            }
        }
        V = V.transpose().eval();
        Eigen::MatrixXcd VExp(3, nHam_sel);
        for (size_t j1=0; j1<3; j1++) {
            VExp(j1, Eigen::all) = V * zedExpF(Eigen::all, j1);
        }
        if (params.incomm) {
            size_t kk = jj % nHklT;
            if (kk < nHklI) {
                K = K1;
            } else if (kk < 2*nHklI) {
                K = K2;
            } else {
                K = cK1;
            }
        }
        if (params.neutron_output) {
            size_t iHklA = jj * 3;
            Eigen::Vector3d hklAN;
            hklAN << inputs.hklA[iHklA], inputs.hklA[iHklA + 1], inputs.hklA[iHklA + 2];
            if (std::isnan(hklAN[0]) || std::isnan(hklAN[1]) || std::isnan(hklAN[2])) {
                if (jj < (params.nHkl - 1)) {
                    hklAN << inputs.hklA[iHklA + 3], inputs.hklA[iHklA + 4], inputs.hklA[iHklA + 5];
                } else {
                    hklAN << 1., 0., 0.;
                }
            }
            qPerp = m1 - (hklAN * hklAN.transpose().eval());
        }
        for (int j1=0; j1<nHam_sel; j1++) {
            Eigen::Map<Eigen::Matrix3cd> Sab(Sab_ptr);
            Sab = VExp(Eigen::all, j1) * VExp(Eigen::all, j1).adjoint();
            if (params.incomm) {
                if (params.helical) {
                    // integrating out the arbitrary initial phase of the helix
                    Eigen::Matrix3cd tmp = (nx * Sab * nx) - ((K2 - m1) * Sab * K2) - (K2 * Sab * (2*K2 - m1));
                    Sab = 0.5 * (Sab - tmp);
                }
                Sab = Sab * K;
            }
            if (params.nTwin > 1) {
                // Rotates the correlation function by the twin rotation matrix
                size_t iT = jj / nHklT;
                Eigen::Map<const Eigen::Matrix3d> rotC(inputs.rotc + iT*9, 3, 3);
                Sab = rotC * Sab * rotC.transpose().eval();
            }
            Sab /= params.nCell;
            if (params.neutron_output) {
                if (params.nformula > 0) {
                    Sab /= params.nformula;
                }
                Sab = (Sab + Sab.transpose().eval()) / 2;
                *(Sperp_ptr++) = (qPerp.cwiseProduct(Sab)).sum();
            } else {
                Sab_ptr += 9;   // This trick only works for column-major data layouts!
            }
        }
    }
    if (params.neutron_output) {
        delete[](Sab_ptr); }
}

template <typename T> T getVal(const mxArray *d) {
    mxClassID cat = mxGetClassID(d);
    switch (cat) {
        case mxDOUBLE_CLASS: return (T)*(mxGetDoubles(d)); break;
        case mxSINGLE_CLASS: return (T)*(mxGetSingles(d)); break;
        case mxINT32_CLASS:  return (T)*(mxGetInt32s(d));  break;
        case mxUINT32_CLASS: return (T)*(mxGetUint32s(d)); break;
        case mxINT64_CLASS:  return (T)*(mxGetInt64s(d));  break;
        case mxUINT64_CLASS: return (T)*(mxGetUint64s(d)); break;
        case mxINT16_CLASS:  return (T)*(mxGetInt16s(d));  break;
        case mxUINT16_CLASS: return (T)*(mxGetUint16s(d)); break;
        case mxINT8_CLASS:   return (T)*(mxGetInt8s(d));   break;
        case mxUINT8_CLASS:  return (T)*(mxGetUint8s(d));  break;
        default:
            throw std::runtime_error("Unknown mxArray class");
    }
}
template <> bool getVal(const mxArray *d) { return (bool)*((mxLogical*)mxGetData(d)); }
template <typename T>
T getField(const mxArray *data_ptr, size_t idx, const char *fieldname, T default_val) {
    mxArray *field = mxGetField(data_ptr, idx, fieldname);
    if (field == nullptr) {
        return default_val;
    }
    return getVal<T>(field);
}

void checkDims(std::string varname, const mxArray *data, size_t dim1, size_t dim2) {
    const mwSize *dims = mxGetDimensions(data);
    if ((size_t)dims[0] != dim1 || (size_t)dims[1] != dim2) {
        throw std::runtime_error("Input " + varname + " has the incorrect size");
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    if (nrhs < 15) {
        throw std::runtime_error("swloop: Requires 15 arguments");
    }
    err_flag.store(0, std::memory_order_release);
    size_t nThreads;
    int nT = getField<int>(prhs[0], 0, "nThreads", -1);
    if (nT < 0) {
        // Gets number of cores / hyperthreads
        nThreads = std::thread::hardware_concurrency() / 2;
        if (nThreads == 0) {
            nThreads = 1;   // Defaults to single thread if not known
        }
    } else {
        nThreads = static_cast<size_t>(nT);
    }
    // Initialises the parameters structure
    struct pars params;
    if (!mxIsStruct(prhs[0])) {
        throw std::runtime_error("swloop: Error first argument must be a param struct");
    }
    params.hermit = getField(prhs[0], 0, "hermit", false);
    params.formfact = getField(prhs[0], 0, "formfact", false);
    params.incomm = getField(prhs[0], 0, "incomm", false);
    params.helical = getField(prhs[0], 0, "helical", false);
    params.bq = getField(prhs[0], 0, "bq", false);
    params.field = getField(prhs[0], 0, "field", false);
    params.fastmode = getField(prhs[0], 0, "fastmode", false);
    params.neutron_output = getField(prhs[0], 0, "neutron_output", false);
    params.omega_tol = getField(prhs[0], 0, "omega_tol", 1.0e-5);
    params.nCell = getField(prhs[0], 0, "nCell", 1.0);
    params.nformula = getField<size_t>(prhs[0], 0, "nformula", 1);
    params.nTwin = getField<size_t>(prhs[0], 0, "nTwin", 1);
    params.nHkl = (size_t)mxGetDimensions(prhs[1])[1];      // hklExt
    params.nBond = (size_t)mxGetDimensions(prhs[5])[1];     // dR
    params.nMagExt = (size_t)mxGetDimensions(prhs[6])[1];   // RR
    params.nBqBond = (size_t)mxGetDimensions(prhs[10])[1];  // bqdR
    // Checks all inputs have the correct dimensions
    size_t len_abcd = 3 * params.nBond;
    size_t nHam = 2 * params.nMagExt;
    checkDims("hklExt", prhs[1], 3, params.nHkl);
    checkDims("ABCD", prhs[2], 1, len_abcd);
    checkDims("idxAll", prhs[3], len_abcd, 2);
    checkDims("ham_diag", prhs[4], nHam, 1);
    checkDims("dR", prhs[5], 3, params.nBond);
    checkDims("RR", prhs[6], 3, params.nMagExt);
    checkDims("S0", prhs[7], 1, params.nMagExt);
    checkDims("zed", prhs[8], 3, params.nMagExt);
    if (params.formfact) {
        checkDims("FF", prhs[9], params.nMagExt, params.nHkl);
    }
    if (params.bq) {
        checkDims("bqdR", prhs[10], 3, params.nBqBond);
        checkDims("bqABCD", prhs[11], 1, 3 * params.nBqBond);
        checkDims("idxBq", prhs[12], 3 * params.nBqBond, 2);
        checkDims("bq_ham_d", prhs[13], nHam, 1);
    }
    // Process other inputs
    struct swinputs inputs;
    if (!mxIsComplex(prhs[8])) {
        throw std::runtime_error("swloop: zed (arg 9) must be complex"); }
    for (size_t ii=1; ii < 15; ii++) {
        if (ii != 2 && ii != 8 && mxIsComplex(prhs[ii])) {
            throw std::runtime_error("swloop: Error an input was found to be complex when it "
                "is expected to be real"); }
    }
    std::complex<double> *ABCDmem;
    if (mxIsComplex(prhs[2])) {
        ABCDmem = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[2]));
    } else {
        ABCDmem = new std::complex<double>[len_abcd];
        double *abcdtmp = mxGetDoubles(prhs[2]);
        for (size_t ii=0; ii<len_abcd; ii++) {
            ABCDmem[ii] = std::complex<double>(abcdtmp[ii], 0.0); }
    }
    inputs.hklExt = mxGetDoubles(prhs[1]);
    inputs.ABCD = ABCDmem;
    inputs.idxAll = mxGetDoubles(prhs[3]);
    inputs.ham_diag = mxGetDoubles(prhs[4]);
    inputs.dR = mxGetDoubles(prhs[5]);
    inputs.RR = mxGetDoubles(prhs[6]);
    inputs.S0 = mxGetDoubles(prhs[7]);
    inputs.zed = reinterpret_cast<cd_t*>(mxGetComplexDoubles(prhs[8]));
    inputs.FF = mxGetDoubles(prhs[9]);
    inputs.bqdR = mxGetDoubles(prhs[10]);
    inputs.bqABCD = mxGetDoubles(prhs[11]);
    inputs.idxBq = mxGetDoubles(prhs[12]);
    inputs.bq_ham_diag = mxGetDoubles(prhs[13]);
    if (params.field) {
        for (size_t ii=0; ii<params.nTwin; ii++) {
            const mxArray *mf_ptr = mxGetCell(prhs[14], ii);
            inputs.ham_MF_v.push_back(mxGetDoubles(mf_ptr));
        }
    }
    inputs.n = mxGetDoubles(mxGetField(prhs[0], 0, "n"));
    inputs.rotc = mxGetDoubles(mxGetField(prhs[0], 0, "rotc"));
    inputs.hklA = mxGetDoubles(mxGetField(prhs[0], 0, "hklA"));

    int *idx0 = new int[params.nHkl];
    double K1[9], K2[9];
    size_t nHklT = params.nHkl / params.nTwin;
    if (params.incomm) {
        // Define index into original Hkl (HklExt is [Q-km, Q, Q+km] - we want origin Q part)
        size_t nHkl0 = nHklT / 3;
        for (size_t jj=0; jj<params.nTwin; jj++) {
            size_t t0 = jj*nHklT;
            for (size_t rr=0; rr<3; rr++) {
                size_t r0 = rr*nHkl0;
                for (size_t kk=0; kk<nHkl0; kk++) {
                    idx0[t0 + kk + r0] = static_cast<int>(t0 + nHkl0 + kk);
                }
            }
        }

    } else {
        for (size_t ii=0; ii<params.nHkl; ii++) {
            idx0[ii] = static_cast<int>(ii); }
    }
    inputs.idx0 = idx0;

    // Creates outputs
    size_t nHam_sel = params.fastmode ? nHam / 2.0 : nHam;
    size_t nn = nHam_sel, nHkl0 = params.nHkl;
    if (params.incomm) {
        nn *= 3;
        nHkl0 = params.nHkl / 3;
    }
    const mwSize dSab[] = {3, 3, nn, nHkl0};
    plhs[0] = mxCreateDoubleMatrix(nn, nHkl0, mxCOMPLEX);
    size_t Sab_sz;
    if (params.neutron_output) {
        plhs[1] = mxCreateDoubleMatrix(nn, nHkl0, mxCOMPLEX);
        Sab_sz = 1;
    } else {
        plhs[1] = mxCreateNumericArray(4, dSab, mxDOUBLE_CLASS, mxCOMPLEX);
        Sab_sz = 9;
    }
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *warn1 = mxGetDoubles(plhs[2]); *warn1 = 0.;
    double *orthwarn = mxGetDoubles(plhs[3]); *orthwarn = 0.;
    struct swoutputs outputs;

    if (params.nHkl < 10 * nThreads || nThreads == 1) {
        // Too few nHkl to run in parallel
        if (params.incomm) {
            outputs.omega = new std::complex<double>[nHam_sel * params.nHkl];
            outputs.Sab = new std::complex<double>[Sab_sz * nHam_sel * params.nHkl];
        } else {
            outputs.omega = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[0]));
            outputs.Sab = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[1]));
        }
        swcalc_fun(0, params.nHkl, std::ref(params), std::ref(inputs), std::ref(outputs));
        int errcode = err_flag.load(std::memory_order_relaxed);
        if (errcode > 0) {
            delete[]idx0;
            if (!mxIsComplex(prhs[2])) {
                delete[]ABCDmem; }
            if (params.incomm) {
                delete [](outputs.omega);
                delete [](outputs.Sab);
            }
            mexErrMsgIdAndTxt("swloop:error", errmsgs[errcode].c_str());
        }
        if (outputs.warn_posdef) *warn1 = 1.;
        if (params.incomm) {
            // Re-arranges incommensurate modes into [Q-km, Q, Q+km] modes
            std::complex<double> *dest_om = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[0]));
            std::complex<double> *dest_Sab = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[1]));
            size_t nHams = nHam_sel * Sab_sz;
            size_t blkSz = nHam_sel * sizeof(std::complex<double>), blkSzs = blkSz * Sab_sz;
            size_t nHklI = params.nHkl / params.nTwin / 3;
            size_t nHklI1 = nHklI * nHam_sel, nHklS1 = nHklI1 * Sab_sz;
            size_t nHklI2 = 2 * nHklI * nHam_sel, nHklS2 = nHklI2 * Sab_sz;
            for (size_t ii=0; ii<params.nTwin; ii++) {
                size_t t0 = ii * nHklT * nHam_sel, s0 = t0 * Sab_sz;
                for (size_t jj=0; jj<nHklI; jj++) {
                    size_t jx = jj * nHam_sel, jy = jx * Sab_sz;
                    memcpy(dest_om, &outputs.omega[jx + t0], blkSz); dest_om += nHam_sel;
                    memcpy(dest_om, &outputs.omega[jx + t0 + nHklI1], blkSz); dest_om += nHam_sel;
                    memcpy(dest_om, &outputs.omega[jx + t0 + nHklI2], blkSz); dest_om += nHam_sel;
                    memcpy(dest_Sab, &outputs.Sab[jy + s0], blkSzs); dest_Sab += nHams;
                    memcpy(dest_Sab, &outputs.Sab[jy + s0 + nHklS1], blkSzs); dest_Sab += nHams;
                    memcpy(dest_Sab, &outputs.Sab[jy + s0 + nHklS2], blkSzs); dest_Sab += nHams;
                }
            }
            delete [](outputs.omega);
            delete [](outputs.Sab);
        }
    } else {
        std::vector<struct swoutputs> outputs_v(nThreads);
        std::vector<std::thread> threads;
        size_t nBlock = params.nHkl / nThreads;
        size_t i0 = 0, i1 = nBlock;
        for (size_t ii=0; ii<nThreads; ii++) {
            outputs_v[ii].omega = new std::complex<double>[nHam_sel * (i1 - i0)];
            outputs_v[ii].Sab = new std::complex<double>[Sab_sz * nHam_sel * (i1 - i0)];
            outputs_v[ii].warn_posdef = false;
            outputs_v[ii].warn_orth = false;
            threads.push_back(
                std::thread(swcalc_fun, i0, i1, std::ref(params), std::ref(inputs), std::ref(outputs_v[ii]))
            );
            i0 = i1;
            i1 += nBlock;
            if (i1 > params.nHkl || ii == (nThreads - 2)) i1 = params.nHkl;
        }
        i0 = 0; i1 = nBlock;
        std::complex<double> *omega_ptr = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[0]));
        std::complex<double> *Sab_ptr = reinterpret_cast<cd_t*>(mxGetComplexDoubles(plhs[1]));

        size_t nHams, blkSz, blkSzs, nHklI;
        if (params.incomm) {
            nHams = nHam_sel * Sab_sz;
            blkSz = nHam_sel * sizeof(std::complex<double>);
            blkSzs = blkSz * Sab_sz;
            nHklI = params.nHkl / params.nTwin / 3;
        }
        for (size_t ii=0; ii<nThreads; ii++) {
            if (threads[ii].joinable()) {
                threads[ii].join();
            }
        }
        int errcode = err_flag.load(std::memory_order_relaxed);
        if (errcode > 0) {
            delete[]idx0;
            if (!mxIsComplex(prhs[2])) {
                delete[]ABCDmem; }
            for (size_t ii=0; ii<nThreads; ii++) {
                delete [](outputs_v[ii].omega);
                delete [](outputs_v[ii].Sab);
            }
            throw std::runtime_error(errmsgs[errcode].c_str());
        }
        for (size_t ii=0; ii<nThreads; ii++) {
            if (params.incomm) {
                // Re-arranges incommensurate modes into [Q-km, Q, Q+km] modes
                for (size_t jj=i0; jj<i1; jj++) {
                    size_t t0 = (jj / nHklT) * nHklT;  // Twin offset
                    size_t kk = (jj % nHklT) / nHklI;  // Flag indicating: 0==Q-km, 1==Q, 2==Q+km
                    size_t k0 = ((jj - t0) - (kk * nHklI))*3 + kk;
                    memcpy(omega_ptr + (t0 + k0) * nHam_sel, &outputs_v[ii].omega[(jj - i0) * nHam_sel], blkSz);
                    memcpy(Sab_ptr + (t0 + k0) * nHams, &outputs_v[ii].Sab[(jj - i0) * nHams], blkSzs);
                }
            } else {
                size_t msz = (i1 - i0) * nHam_sel;
                memcpy(omega_ptr, outputs_v[ii].omega, msz * sizeof(std::complex<double>));
                memcpy(Sab_ptr, outputs_v[ii].Sab, Sab_sz * msz * sizeof(std::complex<double>));
                omega_ptr += msz;
                Sab_ptr += Sab_sz * msz;
            }
            i0 = i1;
            i1 += nBlock;
            if (i1 > params.nHkl || ii == (nThreads - 2)) i1 = params.nHkl;
            if (outputs_v[ii].warn_posdef) *warn1 = 1.;
            delete [](outputs_v[ii].omega);
            delete [](outputs_v[ii].Sab);
        }
    }

    // Clean up
    delete[]idx0;
    if (!mxIsComplex(prhs[2])) {
        delete[]ABCDmem; }
}
