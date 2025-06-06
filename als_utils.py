from time import process_time
from warnings import warn

import numpy as np
from scipy.linalg import svd as sp_svd, solve_discrete_lyapunov, inv as sp_inv, block_diag


def sp_pinv_with_retry(a, atol=None, rtol=None, return_rank=False, check_finite=True, cond=None, rcond=None):
    """
    A reimplemented version of scipy.linalg.pinv that at least tries to get the singular value decomposition through
    gesvd if the default approach through gesdd fails. Gesvd is the default in Octave.
    This is necessary as no other approach, neither using numpy.linalg.pinv (because its svd does not expose the
    option to set the lapack_driver) nor setting the tolerances in exactly the same way as Octave worked to resolve the
    issue that some apparently invertible matrices were not inverted by scipy.linalg.pinv (instead returning a
    LinAlgError) but were inverted by Octave.

    After some searching I found that I am not the only one with this issue, or idea to solve it:
    https://github.com/numpy/numpy/issues/1588
    https://github.com/tenpy/tenpy/blob/main/tenpy/linalg/svd_robust.py
    Those two links also contain some more information on the issue and a better explanation of the problem than I can
    give here.

    Even better: There is a quite old error on Windows that is related to MKL being weird about some random matrices.
    One of them is in tests/float64_pinv_error_matrix.csv. It is a 81x81 matrix with a condition number of 1.45e+16.
    When using scipy's svd implementation or pinv, it will lead to an error occuring on the machine code level with just
    a cryptic error message: "** On entry to DLASCL parameter number  4 had an illegal value". There is no (to me)
    readily apparent feature of that matrix (finiteness, positive definiteness, etc.) that would justify an error on
    that level. The occurrence of the error in this program is also not deterministic, but seems to occur more
    frequently with larger matrices. The error does not occur with Octave's svd implementation.

    Or you can also get "** On entry to SLASCL parameter number 4 had an illegal value" while using float32 arrays.
    This shows that the only way to avoid this is to use OpenBLAS instead of MKL.

    Some links to the error:
    https://stackoverflow.com/a/64738537
    https://stackoverflow.com/a/64971496
    https://github.com/numpy/numpy/issues/16744

    My suspicions:
        - The error is caused by however gesdd is implemented in MKL. Numpy uses gesdd by default, but does not cause
        errors, scipy does too, but I suspect that it is optimized in a different way, since my tests show that it is
        faster than numpy's implementation. Octave uses gesvd, which is slower, but does not cause errors. When using
        Scipy with gesvd, the error does not occur either. Another advantage of gesvd is also its robustness, see above.
        - The error is also disappearing when using the same input, just as a float32 array instead of float64. This
        might be because gesdd is implemented differently for float32 and float64, or because the error is caused by
        some weirdness with the numeric behavior of the more accurate arrays.

    Result: This function is first trying to calculate the pseudoinverse using gesdd, while translating the input array
    to float32. If that fails, it is trying gesvd. The underlying numpy install HAS TO BE compiled with OpenBLAS, not
    MKL, otherwise gesdd will fail with the same error as shown above, either with float32 or float64 arrays.

    See the documentation of scipy.linalg.pinv for more information.
    """
    a = np.asarray(a, dtype=None, order=None)
    if check_finite:
        if a.dtype.char in np.typecodes['AllFloat'] and not np.isfinite(a).all():
            raise ValueError("array must not contain infs or NaNs")
            # added to avoid importing private function _asarray_validated from scipy.linalg._decomp

    try:
        u, s, vh = sp_svd(a, compute_uv=True, full_matrices=False, check_finite=False, lapack_driver='gesdd')
    except np.linalg.LinAlgError:
        u, s, vh = sp_svd(a, compute_uv=True, full_matrices=False, check_finite=False, lapack_driver='gesvd')
        warn('SVD did not converge in gesdd, trying gesvd instead', RuntimeWarning, stacklevel=2)
        # literally all we add... If gesdd fails, we try gesvd. No idea why this is not the default behavior.
        # gesdd is supposed to be faster, but gesvd is the default in Octave and more robust.

    t = u.dtype.char.lower()
    maxS = np.max(s)

    if rcond or cond:
        warn('Use of the "cond" and "rcond" keywords are deprecated and '
             'will be removed in future versions of SciPy. Use "atol" and '
             '"rtol" keywords instead', DeprecationWarning, stacklevel=2)

    # backwards compatible only atol and rtol are both missing
    if (rcond or cond) and (atol is None) and (rtol is None):
        atol = rcond or cond
        rtol = 0.

    atol = 0. if atol is None else atol
    rtol = max(a.shape) * np.finfo(t).eps if (rtol is None) else rtol

    if (atol < 0.) or (rtol < 0.):
        raise ValueError("atol and rtol values must be positive.")

    val = atol + maxS * rtol
    rank = np.sum(s > val)

    u = u[:, :rank]
    u /= s[:rank]
    B = (u @ vh[:rank]).conj().T

    if return_rank:
        return B, rank
    else:
        return B


def autocov_calc(Yi, N):
    """
    Initially generated with SMOP 0.41-beta from inst\autocov_calc.m

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

    Eyy = []
    p, Yi_len = Yi.shape
    for i in range(N):
        temp = Yi[:, i:] @ Yi[:, :Yi.shape[1]-i].T
        temp = temp / (Yi_len - i)
        Eyy.append(temp)
    Eyy = np.vstack(Eyy)
    bhat = Eyy.ravel(order="F")  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    bhat_mtr = np.zeros((N, p ** 2))
    for i in range(p):
        for j in range(p):
            for k in range(N):
                bhat_mtr[k, p * i + j] = Eyy[p * k + i, j]

    npts = Yi_len - N
    cov_bound = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            cov_bound[i, j] = np.sqrt((bhat_mtr[0, p * i + i] * bhat_mtr[0, p * j + j]) / npts)
    cov_bound = cov_bound.ravel(order="F")

    return bhat_mtr, cov_bound, bhat


def autocov_calc_thry(A, C, L, Q, R, N):
    """
    Initially generated with SMOP 0.41-beta from inst\autocov_calc_thry.m

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
    in mind that Octave generally works with column vectors, while numpy uses row vectors by default (for example for the
    output of ravel).
    """

    Abar = A - A @ L @ C
    P = solve_discrete_lyapunov(Abar, Q + A @ L @ R @ L.T @ A.T)
    Eyy = [C @ P @ C.T + R]
    for i in range(N):
        Eyy.append(C @ np.linalg.matrix_power(Abar, i + 1) @ P @ C.T - C @ np.linalg.matrix_power(Abar, i) @ A @ L @ R)
    Eyy = np.vstack(Eyy)
    b = Eyy.ravel(order="F")  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    p = C.shape[0]
    b_mtr = np.zeros((N, p ** 2))
    for i in range(p):
        for j in range(p):
            for k in range(N):
                b_mtr[k, p * i + j] = Eyy[p * k + i, j]

    return b_mtr, b


def comm_mat(m=None, n=None):
    """
    Initially generated with SMOP 0.41-beta from inst\comm_mat.m

    Creates mn x mn commutation matrix Kmn so that for any m x n matrix A,
    Kmn*vec(A) = vec(A')

    Tests run up to m=10 and n=10. Results are identical.
    """
    Kmn = np.zeros((m * n, m * n))

    for i in range(m):
        for j in range(n):
            Kmn[i * n + j, j * m + i] = 1
    return Kmn


def ECM_iter(A, G, C, Q, R):
    """
    Initially generated with SMOP 0.41-beta from inst\ECM_iter.m

    Function to calculate the state error covariance matrix for the DLQE
    using iterations

    Tested. When feeding in the same start parameters the function appears to generate the same output (up to numerical
    precision) as the Octave version.
    """

    p, n = C.shape
    Pold = np.zeros((n, n))
    P = np.eye(n)
    tol = 1e-12
    iter = 0

    while np.linalg.norm(Pold - P, ord='fro') > tol:
        iter += 1
        Pold = P
        P = A @ Pold @ A.T + G @ Q @ G.T - A @ Pold @ C.T @ sp_inv(C @ Pold @ C.T + R) @ C @ Pold @ A.T
        # np.linalg.norm(Pold-P, ord="fro")

    L = P @ C.T @ sp_inv(C @ P @ C.T + R)

    return L, P


def golden_section_Q_mrQ(x00, x33, Q, delQ, A, b, W, nu, vlam_mtr):
    """
    Initially generated with SMOP 0.41-beta from inst\golden_section_Q_mrQ.m

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

    while np.linalg.norm(x3 - x0) > tol:
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


def celldiag(array_iterable, n=0):
    """
    Derived from a function initially generated with SMOP 0.41-beta from inst\celldiag.m

    An alternative implementation of celldiag_old that uses np.pad instead of np.vstack and np.hstack. This is much
    faster for large arrays. This is also not specialized anymore for cell arrays, but works with any iterable of
    arrays.

    Function for diagonalizing matrices

    Tests run, and output compared visually. Results are identical. The original SMOP code appears to return a false
    result with the elements of the cell array in the wrong order. This is not the case for the proper Python
    implementation.

    :param array_iterable: An iterable of arrays
    :type array_iterable: list or tuple
    :param n: Number of rows and columns to pad with zeros. If negative, pad to upper right, if positive, pad to lower
    left.
    :type n: int
    :return: The block diagonal matrix containing elements of array_iterable augmented with zeros if n specified.
    :rtype: np.ndarray
    """
    A = block_diag(*array_iterable)
    n_abs = np.abs(n)
    A = np.pad(A, ((n_abs, 0), (0, n_abs)) if n < 0 else ((0, n_abs), (n_abs, 0)))
    return A


def celldiag_old(F, n=0):
    """
    Initially generated with SMOP 0.41-beta from inst\celldiag.m

    Function for diagonalizing matrices
    function A = celldiag(F,n)
    Function inputs:
    F - cell containing submatrices
    n - number of rows of zeros to add
    if n < 0, zeros are added to upper right
    if n > 0, zeros are added to lower left
    Function outputs:
    A - block diagonal matrix containing elements of F
    (augmented with zeros if n specified)
    A = [0 0; blockdiag(F) 0] if n < 0
    A = [0 blockdiag(F); 0 0] if n > 0
    """
    len_ = np.prod(F.shape)
    A = F[0]
    for i in range(1, len_):
        A = np.vstack([np.pad(A, ((0, 0), (0, F[i].shape[0]))),
                       np.pad(F[i], ((0, 0), (A.shape[0], 0)))])

    n_abs = np.abs(n)
    A = np.pad(A, ((n_abs, 0), (0, n_abs)) if n < 0 else ((0, n_abs), (n_abs, 0)))
    return A


def obj_ls(Q, A, b, W):
    """
    Initially generated with SMOP 0.41-beta from inst\obj_ls.m

    Evaluates only least-squares portion of objective
    (for use in choosing Q from tradeoff curve)

    Tested. The output appears identical to the Octave version up to numerical precision (1e-13).
    """

    y = (A @ Q.ravel(order="F") - b).T @ W @ (A @ Q.ravel(order="F") - b)
    if np.imag(y) != 0:
        y = 10e9
    return y


def obj_tot(Q, A, b, W, nu, vlam_mtr):
    """
    Initially generated with SMOP 0.41-beta from inst\obj_tot.m

    Evaluate objection function value, including penalty on tr(Q) and log-barrier term

    Tests run. Output appears to be identical to the Octave version up to numerical precision (1e-14). I suspect that the
    difference arises in the calculation of the logarithm used in y.
    """
    if np.all(np.linalg.eigvals(Q) > 0):  # Matrix is positive definite, cholesky works
        r = np.linalg.cholesky(Q).T  # Octave calculates the upper Cholesky factor, same as lower (numpy) transposed
        y = (A @ Q.ravel(order="F") - b).T @ W @ (A @ Q.ravel(order="F") - b) + np.trace(
            Q @ vlam_mtr) - nu * 2 * np.sum(np.log(np.diag(r)))
    else:
        y = 1e+100
    return y


def symtran(n):
    """
    Initially generated with SMOP 0.41-beta from inst\symtran.m

    Function to give the duplication matrix that will impose symmetry constraints on
    the covariance matrices, (Q)_s = tran*(Q)_ss

    Tested. The output appears identical to the Octave version.
    """

    r = n * (n + 1) / 2
    tran = np.zeros((n ** 2, int(r)))
    k = 0
    for i in range(n):
        for j in range(i, n):
            temp = np.zeros((n, n))
            temp[i, j] = 1
            if i == j:
                div = 2
            else:
                div = 1
            t2 = (temp + temp.T) / div
            tran[:, k] = t2.ravel(order="F")
            k += 1
    return tran


def symtranT(n=None):
    """
    Initially generated with SMOP 0.41-beta from inst\symtranT.m

    Function to "undo" the symmetry imposing of the duplication matrix
    tran*(x)s = (x)ss
    Equivalent to pinv(symtran(n)), but much faster

    Tested. The output appears identical to the Octave version.
    """

    tran = symtran(n).T
    for i in range(tran.shape[0]):
        tran[i, :] = tran[i, :] / np.sum(tran[i, :])

    return tran


def sdp_QR_mrQ(A, b, W, Q0, lam, nu, pa, Rsym, trstates):
    """
    Initially generated with SMOP 0.41-beta from inst\sdp_QR_mrQ.m

    SemiDefinite programming subroutine.  Uses Newton steps to satisfy
    the KKT conditions with log-barrier function

    Tests run. The difference between this implementation and the Octave implementation appears to be caused by different
    implementations for pinv and cholesky, for example. The numerical precision also appears to play a role.
    Regarding accuracy, see als_sdp_mrQ.
    """

    ntot = Q0.shape[0]
    na = ntot - pa
    ind = np.array(range(1, ntot ** 2 + 1)).reshape((ntot, ntot), order="F")
    # Fortran/column-major, not row-major/C; also 1-based indexing because later we will filter for 0
    _diag_mask = celldiag([np.ones((na, na)), np.ones((pa, pa))])
    ind = (np.logical_and(ind, _diag_mask)) * ind
    Az = np.zeros((A.shape[0], ntot ** 2))

    # Reconstructs least squares matrix for R:
    LH2 = A[:, int(na * (na + 1) / 2):]  # moved here from within if statement below because it is used in both cases
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
    A = np.hstack([A[:, :int(na * (na + 1) / 2)] @ symtranT(na), LH2])
    idx = list(set(ind.ravel(order="F")).difference({0}))
    Az[:, np.subtract(idx, 1)] = A  # return ind to 0-based indexing for Az
    alphmax = 1

    # Initialize
    tol = tol2 = 1e-10
    reltol = 1e-05
    n_iters_tot = 0
    n_iters_max = 100
    Q = Q0
    Qin = sp_inv(Q)
    con2 = 1
    phi = 0
    M = (2 * Az.T) @ W @ Az
    dtemp = (2 * Az.T) @ W @ b

    vlam_mtr = lam * celldiag([np.diag(trstates), np.zeros((pa, pa))])
    vlam = vlam_mtr.ravel(order="F")
    n_iters = 0

    phi_old = 1000.0
    nu_first = 1

    while nu > 1e-12:
        # Decrease penalty on log-barrier term (for sdp constraints)
        if n_iters > 0:
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
            t1 = process_time()
            n_iters_tot += 1
            n_iters += 1
            LH = M + nu * np.kron(Qin, Qin)
            RH = M @ Q.ravel(order="F") - dtemp + vlam - nu * Qin.ravel(order="F")
            delQ = - np.reshape(sp_pinv_with_retry(LH) @ RH, (ntot, ntot), order="F")
            # sp_pinv does not match exactly the implementation in Octave, and also appears to have a different set of
            # tolerances (atol, rtol) instead of octaves single tol. The relative tolerance appears to be calculated
            # by use of the single largest value in x, instead of octaves approach of using norm(x, ord=2) for this.
            # Even when approximating this approach through the setting of atol, rtol sp_pinv returns slightly different
            # results, especially regarding values ~<1e-15, which occur in Octaves output, but not in pythons.
            # This, together with the same happening in the Cholesky decomposition, is the reason for the difference.
            # Update: It appears that this is based on the lapack_driver used. gesvd is the default in Octave, but
            # gesdd is the default in numpy. gesvd is more robust, but slower. gesdd is faster, but less robust.
            # Solved (kind of) by reimplementing sp_pinv_with_retry to use gesvd instead of gesdd if it fails.

            # delQ = -reshape(LH\RH,ntot,ntot); # Can't calculate this way if Hessian is poorly conditioned

            # Optimal Step Size
            alph = golden_section_Q_mrQ(0, alphmax, Q, delQ, Az, b, W, nu, vlam_mtr)
            Qold = Q
            Qnew = Qold + 0.99 * alph * delQ

            if (nu_first < 1) and (min(np.linalg.eig(Qnew)[0]) < 0):  # index 0 is eigenvalues
                Qnew = Qold

            # Update Q
            Q = Qnew
            Q = (Q + Q.T) / 2
            Qin = sp_pinv_with_retry(Q)
            phi_old = phi
            phi = obj_tot(Q, Az, b, W, nu, vlam_mtr)
            con1 = np.abs(phi - phi_old)
            con2 = np.linalg.norm(Q - Qold, ord=2)
            tol1 = tol + reltol * np.abs(phi)
            tol2 = tol + reltol * np.linalg.norm(Q, ord=2)

            t2 = process_time()
            if t2 - t1 > 5:
                # Give warning if code is too slow
                print("WARNING: ALS iterations are really slow for some reason - "
                      "you might want to check what is going on")

        # If solution with the highest penalty on log-barrier term is infeasible, increase log-barrier penalty:
        if nu_first == 1:
            if np.all(np.linalg.eigvals(Q) > 0):  # Matrix is positive definite, cholesky works
                _r = np.linalg.cholesky(Q).T  # Octave calculates the upper Cholesky factor, transpose of the lower
                nu_first = 0
            else:
                # If solution is infeasible, increase nu
                nu *= 100
                n_iters = 0
                Q = Q0
                Qin = sp_inv(Q)
                con2 = 1
                tol2 = tol
                phi = 0

            # Make sure we can't get stuck in this loop
            if nu > 10000000000.0:
                raise RuntimeError("sdp_QR_mrQ: Unable to find >0 solution")

    phi2 = obj_ls(Q, Az, b, W)
    phi_uw = obj_ls(Q, Az, b, np.eye(*W.shape))
    nsq = -nu * np.sum(np.log(np.linalg.eig(Q)[0]))  # index 0 is eigenvalues

    iter_maxed = 0
    if n_iters == n_iters_max:
        iter_maxed = 1

    return Q, phi2, phi, n_iters_tot, nsq, iter_maxed, phi_uw
