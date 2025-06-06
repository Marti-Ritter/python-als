from time import process_time
from warnings import warn

import numba
import numpy as np
from numba import njit, int64, float64, boolean, prange, objmode
from numba.types import Array, Omitted, Tuple
from scipy.linalg import solve_discrete_lyapunov

from .als_utils import sp_pinv_with_retry

arr_f_1d_C = Array(float64, 1, "C")  # Type of numba 1D array of 64-bit floats in C-style (row-major) order
arr_f_2d_C = Array(float64, 2, "C")  # Type of numba 1D array of 64-bit floats in C-style (row-major) order
arr_f_1d_A = Array(float64, 1, "A")  # Type of numba 1D array of 64-bit floats in C-style (row-major) order
arr_f_2d_A = Array(float64, 2, "A")  # Type of numba 1D array of 64-bit floats in C-style (row-major) order
arr_f_2d_F = Array(float64, 2, "F")  # Type of numba 1D array of 64-bit floats in C-style (row-major) order


def first_inv_then_pinv(a):
    try:
        ret = np.linalg.inv(a)
    except np.linalg.LinAlgError:
        warn("Matrix cannot be inverted, using pseudo-inverse instead")
        ret = sp_pinv_with_retry(a, check_finite=True, return_rank=False)
    return ret


@njit([arr_f_2d_C(arr_f_2d_C)], cache=True)
def _numba_pinv(a):
    """
    Since it appears that we can not even use numpy.linalg.pinv as a default as it causes the python kernel to die
    without an error when used in this context in numba we have to go straight for the alternate scipy.linalg.pinv
    implementation. This is a wrapper to make it work with numba. Be mindful that the objmode functionality is
    apparently not stable and barely documented: https://numba.readthedocs.io/en/stable/user/withobjmode.html

    :param a: The array to invert
    :type a: np.ndarray
    :return: The inverted array
    :rtype: np.ndarray
    """
    ret_array = np.zeros_like(a)
    with objmode():
        ret_array[:] = sp_pinv_with_retry(a, check_finite=True, return_rank=False)
    return ret_array


# parallel causes a massive slowdown, so it is disabled
@njit(Tuple((arr_f_2d_C, arr_f_1d_C, arr_f_1d_C))(arr_f_2d_C, int64), parallel=False, cache=True)
def autocov_calc(Yi, N):
    """
    Initially generated with SMOP 0.41-beta from inst/autocov_calc.m

    function [bhat_mtr cov_bound bhat] = autocov_calc(Yi,N)
    Given data Yi, estimates autocovariances up to lag N-1

    Function inputs :
    Yi - Nd x p matrix of data, where Nd is number of samples
    and p is number of measured variables
    N - number of lags

    Function outputs :
    bhat_mtr - N x p matrix of autocovariances
    bhat_mtr =
    [ <y_1(k)y_1(k)>   <y_1(k)y_2(k)>   ...   <y_p-1(k)y_p(k)>     <y_p(k)y_p(k)>
    ...	              ...                     ...                 ...
    <y_1(k)y_1(k-N+1)> <y_1(k)y_2(k-N+1)> ... <y_p-1(k)y_p(k-N+1)> <y_p(k)y_p(k-N+1)>];
    cov_bound - estimated standard deviations of autocovariances
    cov_bound = [sigma_11 sigma_12 ... sigma_1p sigpa_pp]'
    bhat - vector of autocovariances as used in ALS

    Function tested and found to be correct. Returned values are the same up to 1e-15.
    Only difference is that bhat is a column vector in Octave, but a row vector (1-dimensional) in Python.
    bhat_mtr is the same in both cases, a 2-dimensional array with the same values in the same order.
    """

    p, Yi_len = Yi.shape
    Eyy = np.zeros((N*p, p))
    for i in prange(N):
        left_temp = np.ascontiguousarray(Yi[:, i:])
        right_temp = np.ascontiguousarray(Yi[:, :Yi.shape[1]-i].T)
        temp = left_temp @ right_temp
        # temp = Yi[:, i:] @ Yi[:, :Yi.shape[1]-i].T
        temp = temp / (Yi_len - i)
        Eyy[p * i:p * (i + 1), :] = temp
    bhat = Eyy.T.ravel()

    bhat_mtr = np.zeros((N, p ** 2))
    for i in prange(p):
        for j in prange(p):
            for k in prange(N):
                bhat_mtr[k, p * i + j] = Eyy[p * k + i, j]

    npts = Yi_len - N
    cov_bound = np.zeros((p, p))
    for i in prange(p):
        for j in prange(p):
            cov_bound[i, j] = np.sqrt((bhat_mtr[0, p * i + i] * bhat_mtr[0, p * j + j]) / npts)
    cov_bound = cov_bound.T.ravel()

    return bhat_mtr, cov_bound, bhat


@njit(Tuple((arr_f_2d_C, arr_f_1d_C))(arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, int64),
      cache=True)
def _numba_autocov_calc_thry(Abar, P, A, C, L, R, N):
    p = C.shape[0]
    Eyy = np.zeros(((N + 1) * p, p))
    Eyy[:p, :] = C @ P @ C.T + R

    for i in range(N):
        start_idx = (i + 1) * p
        end_idx = (i + 2) * p
        Eyy[start_idx:end_idx, :] = C @ np.linalg.matrix_power(
            Abar, i + 1) @ P @ C.T - C @ np.linalg.matrix_power(Abar, i) @ A @ L @ R
    b = Eyy.T.ravel()  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    b_mtr = np.zeros((N, p ** 2))
    for i in range(p):
        for j in range(p):
            for k in range(N):
                b_mtr[k, p * i + j] = Eyy[p * k + i, j]

    return b_mtr, b


def autocov_calc_thry(A, C, L, Q, R, N):
    """
    Initially generated with SMOP 0.41-beta from inst/autocov_calc_thry.m

    function [b_mtr b] = autocov_calc_thry(A,C,L,Q,R,N)
    Calculates theoretical autocovariances given system matrices
    and noise covariances

    Function inputs :
    System matrices A and C
    Estimator gain L
    Process noise covariance Q (note G is assumed identity)
    Measurement noise covariance R
    Number of lags N

    Function outputs :
    b_mtr - N x p matrix of autocovariances
    b_mtr =
    [ E(y_1(k)y_1(k))   E(y_1(k)y_2(k))   ...   E(y_p-1(k)y_p(k))     E(y_p(k)y_p(k))
        ...	              ...                     ...                 ...
    E(y_1(k)y_1(k-N+1)) E(y_1(k)y_2(k-N+1)) ... E(y_p-1(k)y_p(k-N+1)) E(y_p(k)y_p(k-N+1))];
    cov_bound - estimated standard deviations of autocovariances
    b - vector of autocovariances as used in ALS

    Function tested and found to be correct. Returned values are the same up to 1e-15.
    Only difference is that b is a column vector in Octave, but a row vector (1-dimensional) in Python.

    b_mtr is the same in both cases, a 2-dimensional array with the same values in the same order. For b one has to keep
    in mind that Octave generally works with column vectors, while numpy uses row vectors by default (for example for
    the output of ravel).
    """

    Abar = A - A @ L @ C
    P = solve_discrete_lyapunov(Abar, Q + A @ L @ R @ L.T @ A.T)

    return _numba_autocov_calc_thry(Abar, P, A, C, L, R, N)


# for some reason parallelization leads to a massive increase in time, so it is disabled
@njit(arr_f_2d_C(int64, int64), parallel=False, cache=True)
def comm_mat(m=None, n=None):
    """
    Initially generated with SMOP 0.41-beta from inst/comm_mat.m

    Creates mn x mn commutation matrix Kmn so that for any m x n matrix A,
    Kmn*vec(A) = vec(A')

    Tests run up to m=10 and n=10. Results are identical.
    """
    Kmn = np.zeros((m * n, m * n))

    for i in prange(m):
        for j in prange(n):
            Kmn[i * n + j, j * m + i] = 1
    return Kmn


@njit(Tuple((arr_f_2d_C, arr_f_2d_C))(arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C), cache=True)
def ECM_iter(A, G, C, Q, R):
    """
    Initially generated with SMOP 0.41-beta from inst/ECM_iter.m

    Function to calculate the state error covariance matrix for the DLQE
    using iterations

    Tested. When feeding in the same start parameters the function appears to generate the same output (up to numerical
    precision) as the Octave version.
    """

    p, n = C.shape
    Pold = np.zeros((n, n))
    P = np.eye(n)
    tol = 1e-12

    while np.linalg.norm(Pold - P, ord=None) > tol:  # None is Frobenius norm for matrices (this case) 2-norm for vecs
        Pold = P
        P = A @ Pold @ A.T + G @ Q @ G.T - A @ Pold @ C.T @ np.linalg.inv(C @ Pold @ C.T + R) @ C @ Pold @ A.T

    L = P @ C.T @ np.linalg.inv(C @ P @ C.T + R)

    return L, P


@njit([arr_f_2d_A(arr_f_2d_A, arr_f_2d_A, boolean),
       arr_f_2d_A(arr_f_2d_A, arr_f_2d_A, Omitted(True))], cache=True)
def _celldiag_2_arrays(A, B, diagonal=True):
    x, y = A.shape
    a, b = B.shape

    if diagonal:
        return np.vstack((
            np.hstack((A, np.zeros((x, b)))),
            np.hstack((np.zeros((a, y)), B))
        ))
    else:
        return np.vstack((
            np.hstack((np.zeros((x, b)), A)),
            np.hstack((B, np.zeros((a, y))))
        ))


def celldiag(array_iterable, n=0):
    len_ = len(array_iterable)
    if len_ >= 1:
        A = array_iterable[0]
        for i in range(1, len_):
            A = _celldiag_2_arrays(A, array_iterable[i])
    else:
        A = np.array([[]])

    if n > 0:
        A = _celldiag_2_arrays(A, np.zeros((n, n)), diagonal=False)

    if n < 0:
        n_abs = np.abs(n)
        A = _celldiag_2_arrays(np.zeros((n_abs, n_abs)), A, diagonal=False)
    return A


@njit(float64(arr_f_2d_C, arr_f_2d_C, arr_f_1d_C, arr_f_2d_C), cache=True)
def obj_ls(Q, A, b, W):
    """
    Initially generated with SMOP 0.41-beta from inst/obj_ls.m

    Evaluates only least-squares portion of objective
    (for use in choosing Q from tradeoff curve)

    Tested. The output appears identical to the Octave version up to numerical precision (1e-13).
    """

    y = (A @ Q.T.ravel() - b).T @ W @ (A @ Q.T.ravel() - b)
    if np.imag(y) != 0:
        y = 10e9
    return y


@njit(boolean(arr_f_2d_C), cache=True)
def is_pos_def(A):
    M = A.astype(np.complex128)
    return np.all(np.real(np.linalg.eigvals(M+M.T)) > 0)


@njit(float64(arr_f_2d_C, arr_f_2d_C, arr_f_1d_C, arr_f_2d_C, float64, arr_f_2d_C), cache=True)
def obj_tot(Q, A, b, W, nu, vlam_mtr):
    """
    Initially generated with SMOP 0.41-beta from inst/obj_tot.m

    Evaluate objection function value, including penalty on tr(Q) and log-barrier term

    Tests run. Output appears to be identical to the Octave version up to numerical precision (1e-14). I suspect that the
    difference arises in the calculation of the logarithm used in y.
    """
    if np.any(~np.isfinite(Q)):
        # Cholesky will not work if nonfinite values are present (e.g. NaN, Inf)
        y = 1e+100
        return y

    try:
        # Matrix is positive definite, cholesky works
        r = np.linalg.cholesky(Q).T  # Octave calculates the upper Cholesky factor, same as lower (numpy) transposed
        cholesky_works = True
    except Exception:
        cholesky_works = False
    if cholesky_works:
        y = (A @ Q.T.ravel() - b).T @ W @ (A @ Q.T.ravel() - b) + np.trace(
            Q @ vlam_mtr) - nu * 2 * np.sum(np.log(np.diag(r)))
    else:
        y = 1e+100
    return y


@njit(float64(float64, float64, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_1d_C, arr_f_2d_C, float64, arr_f_2d_C),
      cache=True)
def golden_section_Q_mrQ(x00, x33, Q, delQ, A, b, W, nu, vlam_mtr):
    """
    Initially generated with SMOP 0.41-beta from inst/golden_section_Q_mrQ.m

    Uses golden section method to calculate optimal step size for given search direction
    Objective function: y =  norm(A*Q(:) - b)^2_w + trace(Q*vlam_mtr) - nu*log(det(Q));

    Tested. Results appear identical to the Octave version, at least up to numerical precision.
    """

    tol = 0.001
    x0 = x00
    x3 = x33
    alpha = (3 - np.sqrt(5)) / 2
    x1 = x0 + alpha * (x3 - x0)
    x2 = x3 - alpha * (x3 - x0)

    Q1 = Q + x1 * delQ
    J1 = obj_tot(Q1, A, b, W, nu, vlam_mtr)
    Q2 = Q + x2 * delQ
    J2 = obj_tot(Q2, A, b, W, nu, vlam_mtr)

    while (x3 - x0) > tol:
        # hypothetically this could also be used with vectors and matrices but the way the function is used in the
        # code it is only ever called with scalars
        if J2 < J1:
            x0 = x1
            x1 = x2
            J1 = J2
            x2 = x3 - alpha * (x3 - x0)
            Q2 = Q + x2 * delQ
            J2 = obj_tot(Q2, A, b, W, nu, vlam_mtr)
        else:
            x3 = x2
            x2 = x1
            J2 = J1
            x1 = x0 + alpha * (x3 - x0)
            Q1 = Q + x1 * delQ
            J1 = obj_tot(Q1, A, b, W, nu, vlam_mtr)

    if J1 < J2:
        h = x1
    else:
        h = x2

    return h


@njit(arr_f_2d_C(int64), parallel=True, cache=True)
def symtran(n):
    """
    Initially generated with SMOP 0.41-beta from inst/symtran.m

    Function to give the duplication matrix that will impose symmetry constraints on
    the covariance matrices, (Q)_s = tran*(Q)_ss

    Tested. The output appears identical to the Octave version.
    """

    r = n * (n + 1) / 2
    tran = np.zeros((n ** 2, int(r)))

    for i in prange(n):
        for j in prange(i, n):
            k = (j - i) + (i * (n - i)) + int(i**2/2 + i/2)
            temp = np.zeros((n, n))
            temp[i, j] = 1
            if i == j:
                div = 2
            else:
                div = 1
            t2 = (temp + temp.T) / div
            tran[:, k] = t2.ravel()  # symmetric matrix so no matter which order
    return tran


@njit(arr_f_2d_C(int64), parallel=True, cache=True)
def symtranT(n=None):
    """
    Initially generated with SMOP 0.41-beta from inst/symtranT.m

    Function to "undo" the symmetry imposing of the duplication matrix
    tran*(x)s = (x)ss
    Equivalent to pinv(symtran(n)), but much faster

    Tested. The output appears identical to the Octave version.
    """

    tran = symtran(n).T
    tran_copy = np.zeros(tran.shape)
    for i in prange(tran.shape[0]):
        tran_copy[i, :] = tran[i, :] / np.sum(tran[i, :])
    return tran_copy


@njit
def _typeof(obj):
    with objmode():
        print(obj.shape)
        print(numba.typeof(obj))
    return obj


@njit
def _check_fortran(pos, Qin, Q, Q0, Qnew, delQ):
    with objmode():
        if np.isfortran(Qin):
            print(pos)
            print("Qin")
            Qin = _typeof(Qin)
            print("Q")
            Q = _typeof(Q)
            print("Q0")
            Q0 = _typeof(Q0)
            print("Qnew")
            Qnew = _typeof(Qnew)
            print("delQ")
            delQ = _typeof(delQ)
            print("--------------------")


@njit(float64())
def numba_process_time():
    with objmode(t=float64):
        t = process_time()
    return t


@njit(Tuple((arr_f_2d_C, float64, float64, int64, float64, int64, float64)
            )(arr_f_2d_C, arr_f_1d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, float64, int64, float64),
      cache=True)
def _sdp_QR_mrQ_loop(Az, b, W, Q0, vlam_mtr, nu, ntot, alphmax):
    # Initialize
    vlam = vlam_mtr.ravel()

    tol = tol1 = tol2 = 1e-10
    reltol = 1e-05

    n_iters = 0
    n_iters_tot = 0
    n_iters_max = 100
    iter_maxed = 0

    Q = Q0
    # Numba always returns a fortran array from matrix inversion, so we have to convert it. See here:
    # https://numba.discourse.group/t/numba-and-numpy-have-different-contiguity-results-when-inverting-a-f-order-array/1654
    Qin = np.ascontiguousarray(np.linalg.inv(Q))

    con = con1 = con2 = 1
    phi = 0
    M = (2 * Az.T) @ W @ Az  # Hessian part of norm fit
    dtemp = (2 * Az.T) @ W @ b  # data to be used for gradient

    phi_old = 1e3
    nu_first = True

    try:
        while nu > 1e-12:
            if n_iters > 0:  # Decrease penalty on log-barrier term (for sdp constraints)
                if n_iters / n_iters_max < 0.75:
                    nu *= 0.75
                else:
                    nu *= 0.95 * n_iters / n_iters_max
            n_iters = 1

            phi = obj_tot(Q, Az, b, W, nu, vlam_mtr)  # Objective function value
            con1 = np.abs(phi - phi_old)  # Stopping criterion 1
            tol1 = tol + reltol * np.abs(phi)  # Tolerance for stopping criterion 1

            # Optimize for each penalty on log-barrier term
            while ((con1 > tol1) or (con2 > tol2)) and (n_iters < n_iters_max):
                # see comment below for why this is commented out
                # t1 = numba_process_time()
                n_iters_tot += 1
                n_iters += 1

                # LH = M + nu * np.kron(Qin, Qin)
                # RH = M @ Q.ravel() - dtemp + vlam - nu * Qin.ravel()
                delQ = - (_numba_pinv(M + nu * np.kron(Qin, Qin)) @ (
                            M @ Q.ravel() - dtemp + vlam - nu * Qin.ravel())).reshape((ntot, ntot))

                # Optimal Step Size
                alph = golden_section_Q_mrQ(0, alphmax, Q, delQ, Az, b, W, nu, vlam_mtr)
                Qold = Q
                Qnew = Qold + 0.99 * alph * delQ  # empirical tests show that Qnew is (nearly?) always real symmetric

                if not nu_first and np.min(np.linalg.eigvalsh(Qnew)) < 0:
                    # if not nu_first and np.min(np.real(np.linalg.eigvals(Qnew.astype(np.complex128)))) < 0:
                    Qnew = Qold

                # Update Q
                Q = Qnew
                Q = (Q + Q.T) / 2  # Q must be symmetric!
                Qin = _numba_pinv(Q)  # Qin, as the inverse of Q, must also be symmetric
                phi_old = phi
                phi = obj_tot(Q, Az, b, W, nu, vlam_mtr)
                con1 = np.abs(phi - phi_old)
                con2 = np.linalg.norm(Q - Qold, ord=2)
                tol1 = tol + reltol * np.abs(phi)
                tol2 = tol + reltol * np.linalg.norm(Q, ord=2)

                # if the numba code is ever too slow, we have other problems
                # There is no proper way to implement this without losing the compilation speedup for the outer loop
                """
                t2 = numba_process_time()
                if t2 - t1 > 5:
                    # Give warning if code is too slow
                    print("WARNING: ALS iterations are really slow for some reason - "
                          "you might want to check what is going on")
                """
            # If solution with the highest penalty on log-barrier term is infeasible, increase log-barrier penalty:
            if nu_first:
                try:
                    # Matrix is positive definite, cholesky works
                    _r = np.linalg.cholesky(Q).T  # Octave calculates the upper Cholesky factor, transpose of the lower
                    cholesky_works = True
                except Exception:
                    # excepting every Exception is all numba can handle, but the above try is literally just
                    # a cholesky, so this can only be a LinalgError anyway
                    cholesky_works = False

                if cholesky_works:
                    nu_first = False
                else:
                    # If solution is infeasible, increase nu
                    nu *= 100
                    n_iters = 0
                    Q = Q0
                    Qin = np.ascontiguousarray(np.linalg.inv(Q))
                    con2 = 1
                    tol2 = tol
                    phi = 0

                # Make sure we can't get stuck in this loop
                if nu > 1e10:
                    print("sdp_QR_mrQ: Unable to find >0 solution")
                    raise Exception

    except Exception:
        print("An unexpected error occurred in sdp_QR_mrQ. "
              "The error was caught and the function will return the current values. "
              "Please check your input parameters. If the previous message indicated that no solution was found, "
              "you may want to try different initial values for Q0.")
        iter_maxed = -1  # flag to indicate that an unexpected error occurred

    phi2 = obj_ls(Q, Az, b, W)
    phi_uw = obj_ls(Q, Az, b, np.eye(*W.shape))
    nsq = -nu * np.sum(np.log(np.real(np.linalg.eig(Q.astype(np.complex128))[0])))  # index 0 is eigenvalues

    if (iter_maxed != -1) and (n_iters == n_iters_max):  # check if maximum iters reached and no exception occurred
        # the loop was exited because of the iteration limit
        iter_maxed = 1

    return Q, phi2, phi, n_iters_tot, nsq, iter_maxed, phi_uw


@njit(Tuple((arr_f_2d_C, float64, float64, int64, float64, int64, float64)
            )(arr_f_2d_C, arr_f_1d_C, arr_f_2d_C, arr_f_2d_C, float64, float64, int64, int64, arr_f_1d_C), cache=True)
def sdp_QR_mrQ(A, b, W, Q0, lam, nu, pa, Rsym, trstates):
    """
    Initially generated with SMOP 0.41-beta from inst/sdp_QR_mrQ.m

    SemiDefinite programming subroutine.  Uses Newton steps to satisfy
    the KKT conditions with log-barrier function

    Tests run. The difference between this implementation and the Octave implementation appears to be caused by different
    implementations for pinv and cholesky, for example. The numerical precision also appears to play a role.
    Regarding accuracy, see als_sdp_mrQ.
    """

    ntot = Q0.shape[0]
    na = ntot - pa
    ind = np.arange(1, ntot ** 2 + 1).reshape((ntot, ntot))
    # ind = np.array(range(1, ntot ** 2 + 1)).reshape((ntot, ntot))
    # Fortran/column-major, not row-major/C; also 1-based indexing because later we will filter for 0
    _diag_mask = _celldiag_2_arrays(np.ones((na, na)), np.ones((pa, pa)), diagonal=True)
    # _diag_mask = celldiag([np.ones((na, na)), np.ones((pa, pa))])
    ind = (np.logical_and(ind, _diag_mask)) * ind
    Az = np.zeros((A.shape[0], ntot ** 2))

    # Reconstructs least squares matrix for R:

    # moved here from within if statement below because it is used in both cases
    LH2 = np.ascontiguousarray(A[:, int(na * (na + 1) / 2):])
    
    if pa != 0:
        if Rsym == 1:
            # For symmetric R matrix
            mtrx = symtranT(pa)
        else:
            # For diagonal R matrix
            mtrx = np.zeros((pa, pa ** 2))
            for i in range(pa):
                mtrx[i, i * (pa + 1)] = 1
        LH2 = LH2 @ mtrx

    # Reconstructs least-squares matrix for Q:
    A_slice_contig = np.ascontiguousarray(A[:, :int(na * (na + 1) / 2)])
    A = np.hstack((A_slice_contig @ symtranT(na), LH2))
    idx = ind.ravel()
    idx = idx[idx != 0]
    Az[:, np.subtract(idx, 1)] = A  # return ind to 0-based indexing for Az
    alphmax = 1

    vlam_mtr = lam * _celldiag_2_arrays(np.diag(trstates), np.zeros((pa, pa)), diagonal=True)
    # vlam_mtr = lam * celldiag([np.diag(trstates), np.zeros((pa, pa))])

    # loop extracted to separate function for numba compatibility
    return _sdp_QR_mrQ_loop(Az, b, W, Q0, vlam_mtr, nu, ntot, alphmax)


@njit([arr_f_2d_C(int64, int64, arr_f_1d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_F, arr_f_2d_F),
       arr_f_2d_C(int64, int64, arr_f_1d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_C, arr_f_2d_F)],
      cache=True)
def _estimate_states(na, datapts, xhat0, Aa, Ba, Ca, L, y, u):
    """
    This looks like this because all matrices should be C-contiguous but the original implementation was writing columns
    of xhat_ to the output matrix, which is Fortran-ordered. So we transpose everything, do the calculation, and
    transpose back.
    This is not much, but at least it is a bit faster than the original implementation.
    """
    xhat = np.zeros((datapts, na)).T
    xhat_ = np.zeros((datapts + 1, na)).T
    xhat_[:na, 0] = xhat0
    for i in prange(datapts):
        xhat[:, i] = xhat_[:, i] + L @ (y[:, i] - Ca @ xhat_[:, i])
        xhat_[:, i + 1] = Aa @ xhat[:, i] + Ba @ u[:, i]
    xhat_ = np.ascontiguousarray(xhat_[:, :datapts])
    return xhat_


@njit([Tuple((arr_f_1d_C, arr_f_2d_C))(int64, int64, arr_f_2d_C, int64, int64, Omitted(1)),
       Tuple((arr_f_1d_C, arr_f_2d_C))(int64, int64, arr_f_2d_C, int64, int64, int64)], parallel=True, cache=True)
def _calculate_autocorrelations(datapts, start, inntrun, N, pa, dataweight=1):
    datatrun = datapts - start
    inntrun_width = inntrun.shape[1]

    Eyy = np.zeros((N * pa, pa))
    # Eyy = []
    for i in prange(N):
        _ltemp = np.ascontiguousarray(inntrun[:, i:])
        _rtemp = np.ascontiguousarray(inntrun[:, :inntrun_width - i].T)
        temp = _ltemp @ _rtemp
        temp = temp / (datatrun - i)
        Eyy[pa * i:pa * (i + 1), :] = temp
        # Eyy.append(temp)
    # Eyy = np.vstack(Eyy)
    Eyy = Eyy.T.ravel()  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    if dataweight == 1:
        # Calculation of Weighting matrix for one column ALS
        Nd = inntrun.shape[1]
        nt = Nd - N + 1
        trials = 2 * N
        datapts2 = nt // trials
        covb = np.zeros((N * pa ** 2, N * pa ** 2))

        for i_trial in prange(trials):
            Yysmall = np.zeros((N * pa, datapts2))
            for i in prange(datapts2):
                yyst = inntrun[:, i * trials + i_trial:i * trials + i_trial + N]
                Yysmall[:, i] = yyst.T.ravel()  # see above, octave uses Fortran

            Py = np.cov(Yysmall.T, rowvar=False)  # octave default is column
            Px = np.ascontiguousarray(Py[:pa, :pa])
            Pyx = np.ascontiguousarray(Py[:, :pa])
            covb += 1 / datapts2 * (np.kron(Px, Py) + comm_mat(pa, N * pa) @ np.kron(Pyx, Pyx.T))

        covb /= trials
        Wm = _numba_pinv(covb)
    else:
        # Use identity-based weighting
        Wm = np.eye(len(Eyy))
    return Eyy, np.ascontiguousarray(Wm)


@njit(Tuple((arr_f_2d_C, arr_f_2d_C, arr_f_1d_C))(int64, arr_f_2d_C, int64, arr_f_2d_C, arr_f_2d_C, arr_f_1d_C, int64,
                                                  int64, int64, int64), parallel=True, cache=True)
def _form_autocovariance_matrices(Rsym, LHSsingc, ga, Qest_cell0, Rest_cell0, Eyy, N, pa, datapts, start):
    if Rsym == 1:
        Eyy_thry = LHSsingc @ np.hstack((symtranT(ga) @ Qest_cell0.T.ravel(),
                                         symtranT(pa) @ Rest_cell0.T.ravel()))
    else:
        Eyy_thry = LHSsingc @ np.hstack((symtranT(ga) @ Qest_cell0.T.ravel(), np.diag(Rest_cell0)))

    Eyy = Eyy.T.reshape(pa, N * pa).T
    Eyy_thry = Eyy_thry.T.reshape(pa, N * pa).T

    bhat_mtr = np.zeros((N, pa ** 2))
    bhat_mtr_thry = np.zeros((N, pa ** 2))
    for i in prange(pa):
        for j in prange(pa):
            for k in prange(N):
                bhat_mtr[k, pa * i + j] = Eyy[pa * k + i, j]
                bhat_mtr_thry[k, pa * i + j] = Eyy_thry[pa * k + i, j]

    npts = datapts - start - N
    cov_bound = np.zeros((pa, pa))
    for i in prange(pa):
        for j in prange(pa):
            cov_bound[i, j] = np.sqrt((bhat_mtr[0, pa * i + i] * bhat_mtr[0, pa * j + j]) / npts)
    cov_bound = cov_bound.T.ravel()
    return bhat_mtr, bhat_mtr_thry, cov_bound
