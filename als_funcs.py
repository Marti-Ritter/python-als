from time import process_time

import numpy as np
from control import dlqe
from matplotlib import pyplot as plt
from qpsolvers import solve_qp as qp
from scipy.linalg import solve_discrete_lyapunov

from .als_utils import symtran, ECM_iter, comm_mat, sp_pinv_with_retry, celldiag, symtranT, sdp_QR_mrQ


def als_condno(N, model, estimator, rform=None):
    """
    Initially generated with SMOP 0.41-beta from inst\als_condno.m

    Returns condition number and largest and smaller singular values for
    the ALS matrix (script A) formed from N, model, estimator.  Does not
    solve ALS problem
    
    [cond_no svals Mcond] = als_condno(N,model,estimator,Rform)
    
    Function Inputs :
    N (window size),
    model.A,
    model.C,
    model.G (optional, default is identity matrix)
    estimator.L (initial estimator gain - optional)
    or estimator.Q and estimator.R (process and measurement noise
    covariances, used to calculate estimator gain),
    Rform (optional input to specify structure of R_v.  If
    Rform = "sym", R_v is symmetric, if Rform = "diag", R_v is
    diagonal.  Default is diagonal)
    
    Function outputs :
    cond_no (condition number of script A)
    svals (singular values of script A)
    Mcond (condition number of M = kron(C,I)*(I-kron(A,A))^-1*kron(G,G)*D_g

    For easily invertible matrices (low condition number), the results are the same, but for matrices with a high condition
    number, the results differ. Both implementations claim to use the 2-norm, but the Octave implementation returns a
    different result than the Python implementation.
    Up to the calculation of M both implementations contain identical values, up to numerical precision.
    For the singular value vector svals, the values are the same up to numerical precision (~1e-15).
    """
    Aa = model.A
    Ca = model.C

    if hasattr(model, 'G'):
        Ga = model.G
    else:
        Ga = np.eye(*Aa.shape)

    n, g = Ga.shape
    p = Ca.shape[0]

    if hasattr(estimator, 'L'):
        L = estimator.L
    else:
        L, _P, _E = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R)

    na = n
    ga = g
    pa = p

    Rsym = 0
    if rform is not None:
        if rform == 'sym':
            Rsym = 1
        elif rform == 'diag':
            Rsym = 0
        else:
            print('WARNING: als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal')

    ##################################################
    # Building the constant matrix for the LS problem
    ##################################################

    Ain = Aa - Aa @ L @ Ca
    OO = []
    temp = np.eye(n)

    for i in range(N):
        OO.append(Ca @ temp)
        temp = temp @ Ain
    OO = np.vstack(OO)

    M1 = np.zeros((na ** 2, ga ** 2))
    i = 0
    for j in range(ga):
        for k in range(ga):
            II = np.zeros((ga, ga))
            II[k, j] = 1
            t1 = solve_discrete_lyapunov(Ain, Ga @ II @ Ga.T)
            M1[:, i] = t1.ravel(order="F")
            i += 1

    M2 = np.zeros((na ** 2, pa ** 2))
    i = 0
    for j in range(pa):
        for k in range(pa):
            II = np.zeros((pa, pa))
            II[k, j] = 1
            t2 = solve_discrete_lyapunov(Ain, Aa @ L @ II @ L.T @ Aa.T)
            M2[:, i] = t2.ravel(order="F")
            i += 1

    ##########################################
    # Single column ALS method
    ##########################################

    PSI = np.eye(pa)
    for i in range(N - 1):
        PSI = np.vstack([PSI, -Ca @ np.linalg.matrix_power(Ain, i) @ Aa @ L])

    OOtemp = np.kron(Ca, OO)
    PSItemp = np.kron(np.eye(pa), PSI)

    LH1 = OOtemp @ M1
    LH2 = OOtemp @ M2 + PSItemp

    if Rsym == 1:
        LHSsingc = np.hstack([LH1 @ symtran(ga), LH2 @ symtran(pa)])
    else:
        LHSsingc = np.hstack([LH1 @ symtran(ga), LH2[:, :pa ** 2:pa + 1]])

    M = np.kron(Ca, np.eye(n)) @ M1 @ symtran(ga)

    Mcond = np.linalg.cond(M, p=2)
    cond_no = np.linalg.cond(LHSsingc)
    svals = np.linalg.svd(LHSsingc)[1]

    return cond_no, svals, Mcond


def als_diag(data, N, model, estimator):
    """
    Initially generated with SMOP 0.41-beta from inst\als_diag.m

    The Diagonal ALS function.  Uses Diagonal ALS form to estimate only
    the diagonal elements of Qw and Rv.
    Function Inputs : data.yk (measurements), data.uk (inputs),
    data.xhatk (state estimates- optional), data.datapts (number of data
    points considered), data.start (data to be ignored in the beginning
    till initial condition is negligible),
    model.A, model.B (optional), model.C, model.G , N (window size),
    estimator.Q (Qw initial guess), estimator.R (Rv initial guess),
    estimator.L (initial estimator gain - optional).
    Function outputs : estimated Qw, estimated Rv, estimated filter gain
    L, and ALS LHS matrix As and RHS vector bhat.
    Please see file simulate_data8_diag.m for an example of the
    Diagonal ALS implementation for a simple system.
    For Matlab implementation, the call of qp must be replaced by
    quadprog - see line 133


    Tests run. Results show slight deviation between Python and Octave implementations, which is likely due to the different
    quadratic programming solvers used. Its input arguments are the same as those of the Octave function, but it does not
    result in exactly the same output.
    The results are still very close, so the Python implementation is considered to be correct.
    """

    datapts = data.datapts
    Aa = model.A
    Ca = model.C

    if hasattr(model, 'G'):
        Ga = model.G
    else:
        Ga = np.eye(*Aa.shape)

    n, g = Ga.shape
    p = Ca.shape[0]

    if hasattr(model, 'B'):
        Ba = model.B
    else:
        Ba = np.zeros((n, n))
        data.uk = np.zeros((n, datapts))

    if hasattr(data, 'start'):
        start = data.start
    else:
        start = 100

    na = n
    ga = g
    pa = p

    # Estimator Simulation
    y = data.yk
    u = data.uk
    if hasattr(estimator, 'L'):
        L = estimator.L
    else:
        L, _P, _E = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R, method="scipy")

    if not hasattr(data, 'xhatk'):
        xhat = np.zeros((na, datapts))
        xhat_ = np.zeros((n, datapts + 1))
        xhat_[:n, 0] = model.xhat0
        for i in range(datapts):
            xhat[:, i] = xhat_[:, i] + L @ (y[:, i] - Ca @ xhat_[:, i])
            xhat_[:, i + 1] = Aa @ xhat[:, i] + Ba @ u[:, i]
        xhat_ = xhat_[:, :datapts]
    else:
        xhat_ = data.xhatk

    # Estimate L-innovations
    inntrun = y[:, start:] - Ca @ xhat_[:, start:]
    datatrun = datapts - start

    # Calculation of Autocorrelations for one column ALS
    Eyy = []
    inntrun_width = inntrun.shape[1]

    for i in range(N):
        temp = inntrun[:, i:] @ inntrun[:, :inntrun_width - i].T
        temp = temp / (datatrun - i)
        Eyy.append(temp)
    Eyy = np.vstack(Eyy)
    Eyy = Eyy.ravel(order="F")  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    ##################################################
    # Building the constant matrix for the LS problem
    ##################################################

    Ain = Aa - Aa @ L @ Ca
    OO = []
    temp = np.eye(n)
    for i in range(N):
        OO.append(Ca @ temp)
        temp = temp @ Ain
    OO = np.vstack(OO)

    # temporary variables
    M1 = np.zeros((na ** 2, ga ** 2))
    i = 0
    for j in range(ga):
        for k in range(ga):
            II = np.zeros((ga, ga))
            II[k, j] = 1
            t1 = solve_discrete_lyapunov(Ain, Ga @ II @ Ga.T)
            M1[:, i] = t1.ravel(order="F")
            i += 1

    M2 = np.zeros((na ** 2, pa ** 2))
    i = 0
    for j in range(pa):
        for k in range(pa):
            II = np.zeros((pa, pa))
            II[k, j] = 1
            t2 = solve_discrete_lyapunov(Ain, Aa @ L @ II @ L.T @ Aa.T)
            M2[:, i] = t2.ravel(order="F")
            i += 1

    #############################
    # Diagonal ALS method
    #############################

    PSI = np.eye(pa)
    for i in range(N - 1):
        PSI = np.vstack([PSI, -Ca @ np.linalg.matrix_power(Ain, i) @ Aa @ L])

    OOtemp = np.kron(Ca, OO)
    PSItemp = np.kron(np.eye(pa), PSI)

    LH1 = OOtemp @ M1
    LH2 = OOtemp @ M2 + PSItemp

    As_diag = np.hstack([LH1[:, :ga ** 2:ga + 1], LH2[:, :pa ** 2:pa + 1]])

    # Testing the uniqueness of covariance estimates
    Arank = np.linalg.matrix_rank(As_diag, 0.0001)
    nr, nc = As_diag.shape

    if nc > Arank:
        print("WARNING: Covariance estimates are not unique!")

    # original arguments:
    oct_H = As_diag.T @ As_diag
    oct_q = -As_diag.T @ Eyy
    oct_lb = np.zeros(ga + pa)

    # oct_Ain = np.eye(pa + ga)
    # oct_x0 = np.ones(ga + pa)
    # oct_A = oct_b = oct_ub = oct_Alb = oct_Aub = []

    # solving general formula: min 0.5 x'*H*x + x'*q
    # subject to: A*x <= b, lb <= x <= ub, A_lb <= A_in*x <= A_ub

    # qpsolvers formula: min 0.5 x'*P*x + q'*x
    # subject to: G*x <= h, A*x = b, lb <= x <= ub

    # leading to the following substitutions:
    # inequality: G = A_in, h = b; equality: A = A, b = b; vector bounds: lb = lb, ub = ub;
    # matrix bounds: [] = A_lb, [] = A_ub (have no equivalent in qpsolvers); initvals = x0

    # In effect:
    # cost matrix P = H, cost vector q = q, lower bounds lb = lb
    # for some unknowable reason, the HiGHS algorithm returns nearly the same output, even without feeding in the
    # inequality matrix G (oct_Ain, which is the identity matrix in this case) and the inequality vector h (which is not
    # given as an argument to qp in Octave). It also does not take initial values as it does not permit a "warm start"
    # (which would be the only reason to feed in initial values in the first place), so initvals=oct_x0 is left out.

    Xest_diag = qp(P=oct_H, q=oct_q, lb=oct_lb, solver="highs")

    if np.prod(Xest_diag) == 0:
        print("WARNING: Covariance estimate(s) is (are) at constraints! You may have bad data!\n")

    Qest = np.diag(Xest_diag[:ga])
    Rest = np.diag(Xest_diag[ga:])
    Lest = ECM_iter(Aa, Ga, Ca, Qest, Rest)[0]
    As = As_diag
    bhat = Eyy

    return Qest, Rest, Lest, As, bhat


def als_sdp_mrQ(data, N, model, estimator, rho_values=None, tracestates=None, rform=None,
                weight=None, plot=True, use_sdp=True):
    """
    Initially generated with SMOP 0.41-beta from inst\als_sdp_mrQ.m

    A modified ALS-SDP function.  Uses Single Column ALS form and imposes
    semidefinite constraints.
    Options allow semidefinite constraints to be removed.

    [Qest_cell,Rest_cell,trQ,Phi,Phi_tot,phi0,bhat_mtr,bhat_mtr_thry,cov_bound, ...
    Iter,Iter_maxed,timespent] = als_sdp_mrQ(data,N,model,estimator,varargin)

    Function Inputs:

    data.yk (measurements),
    data.uk (inputs),
    data.xhatk (state estimates- optional),
    data.datapts (number of data points considered, optional),
    data.start (data to be ignored in the beginning until initial condition is
    negligible - optional, default is 100),
    N (window size),
    model.A,
    model.B (optional),
    model.C,
    model.G (optional, default is identity matrix)
    estimator.Q (Qw initial guess, optional)
    estimator.R (Rv initial guess, optional),
    estimator.L (initial estimator gain - optional)

    Either (Q,R) or L must be supplied.  If L is supplied, (Q,R) are used as
    initial guess but not in calculating the estimator gain.  If L is not supplied
    then (Q,R) are used to calculate initial estimator

    Optional inputs for vargin:
    'rho_values', [rho1 rho2 rho3 ...] (vector of penalties on trace(Q), default is
    [0])
    'Rform', 'sym' or 'Rform', 'diag' (imposes symmetric or diagonal constraints on
    Rv, default is diagonal)
    'Weight', 'I', or "Weight', 'data' (imposes identity matrix weighting or
    data-based weighting on least-squares term, default is data-based)
    'Plot', 0 (turns off all plots) or 'plot', 1 (produces plots of L-innovations
    and plots of autocovariance estimates and fits), default is to produce plots
    'sdp', 0 (solves unconstrained problem), default is 'sdp', 1 (enforces
    semidefinite constraints)
    'tracestates', [1 1 ... 0 0] (vector indicating which states are to be
    penalized in the trace(Q) term, default is [1 1 ... 1], to penalize all states)


    Function outputs :

    Qest_cell - cell containing estimated Qw for each penalty on tr(Q)
    Rest_cell - cell containing estimated Rv for each penalty on tr(Q)
    trQ - vector containing trace(Q) for each penalty on tr(Q)
    Phi - vector containing least-squares portion of objective for each
    penalty on tr(Q); scaled so that phi = 1 when rho = 0
    Phi_tot - vector containing full objective function value for each
    penalty on tr(Q)
    phi0 - vector containing objective function value at rho = 0 (prior
    to scaling)
    bhat_mtr - matrix of autocovariances estimated from data, columns are
    <y_1(k),y_1(k-j)>, ..., <y_1(k),y_p(k-j)>, ..., <y_p(k),y_1(k-j)>, ...,
    <y_p(k),y_p(k-j)>
    bhat_mtr - matrix of theoretical autocovariances calculated from model
    and ALS results (for smallest penalty on trace(Q))
    cov_bound - estimated standard deviations for each auto-cross covariances,
    [sigma_11; sigma_12; ... sigma_pp], use +/- 2*sigma as 95# confidence intervals
    Iter - vector containing number of iterations for each penalty on tr(Q)
    Iter_maxed(i) = 1 if maximum number of iterations are reached for that
    penalty on trQ
    Timespent - vector of time spent in seconds of CPU time for each
    penalty on tr(Q)

    Tests run. With an example dataset with a rho logspace -6 to 6, 25 values, there were no distinct differences beyond
    the numerically reasonable visible.
    It is my impression that this is accurate enough.
    """
    if not hasattr(data, 'datapts'):
        data.datapts = data.yk.shape[1]

    datapts = data.datapts
    Aa = model.A
    Ca = model.C

    if hasattr(model, 'G'):
        Ga = model.G
    else:
        Ga = np.eye(*Aa.shape)

    n, g = Ga.shape
    p = Ca.shape[0]
    if hasattr(model, 'B'):
        Ba = model.B
    else:
        Ba = np.zeros((n, n))
        data.uk = np.zeros((n, datapts))

    if hasattr(data, 'start'):
        start = data.start
    else:
        start = 100

    na = n
    ga = g
    pa = p

    # Deal with variable arguments
    lam_vec = rho_values
    trstates = tracestates if tracestates is not None else np.ones((ga, 1))
    trstates = np.atleast_1d(trstates.squeeze())  # to ensure it is a row vector, or a matrix

    Rsym = 0
    if rform is not None:
        if rform == 'sym':
            Rsym = 1
        elif rform == 'diag':
            Rsym = 0
        else:
            print('WARNING: als_sdp_mrQ: Unknown structure type for R; defaulting to diagonal')

    dataweight = 1
    if weight is not None:
        if weight == 'I':
            dataweight = 0
        elif weight == 'data':
            dataweight = 1
        else:
            print('WARNING: als_sdp_mrQ: Unknown weighting type; defaulting to data-based')

    plot_flag = plot
    sdp = use_sdp

    # Estimator Simulation
    y = data.yk
    u = data.uk

    # Calculate estimator gain if not supplied:
    if hasattr(estimator, 'L'):
        L = estimator.L
    else:
        L, _P, _E = dlqe(Aa, Ga, Ca, estimator.Q, estimator.R, method="scipy")

    # Estimate states
    if not hasattr(data, 'xhatk'):
        xhat = np.zeros((na, datapts))
        xhat_ = np.zeros((n, datapts + 1))
        xhat_[:n, 0] = model.xhat0
        for i in range(datapts):
            xhat[:, [i]] = xhat_[:, [i]] + L @ (y[:, [i]] - Ca @ xhat_[:, [i]])
            xhat_[:, [i + 1]] = Aa @ xhat[:, [i]] + Ba @ u[:, [i]]
        xhat_ = xhat_[:, :datapts]
    else:
        xhat_ = data.xhatk

    # Estimate L-innovations
    inntrun = y[:, start:] - Ca @ xhat_[:, start:]

    # Optional - use to look at plots to check if error looks white:
    if plot_flag:
        if p > 1:
            n_cols = 2
            n_rows = p // 2

            if p % 2 == 1:
                n_rows += 1
        else:
            n_cols = 1
            n_rows = 1

        fig = plt.figure()

        for i in range(p):
            ax = fig.add_subplot(n_rows, n_cols, i + 1)
            ax.plot(inntrun[i, :], 'm', linewidth=2)
            ax.set_xlabel('Time')
            ax.set_ylabel('Innovation of y_{}'.format(i))
        plt.tight_layout()
        plt.show()

    # Calculation of Autocorrelations for one column ALS
    datatrun = datapts - start
    Eyy = []
    inntrun_width = inntrun.shape[1]

    for i in range(N):
        temp = inntrun[:, i:] @ inntrun[:, :inntrun_width - i].T
        temp = temp / (datatrun - i)
        Eyy.append(temp)
    Eyy = np.vstack(Eyy)
    Eyy = Eyy.ravel(order="F")  # Octave uses Fortran (column-major) order, numpy by default uses C (row-major)

    if dataweight == 1:
        # Calculation of Weighting matrix for one column ALS
        Nd = inntrun.shape[1]
        nt = Nd - N + 1
        trials = 2 * N
        datapts2 = nt // trials
        covb = np.zeros((N * pa ** 2, N * pa ** 2))

        for i_trial in range(trials):
            Yysmall = np.zeros((N * p, datapts2))
            for i in range(datapts2):
                yyst = inntrun[:, i * trials + i_trial:i * trials + i_trial + N]
                Yysmall[:, i] = yyst.ravel(order="F")  # see above, octave uses Fortran

            Py = np.cov(Yysmall.T, rowvar=False)  # octave default is column
            Px = Py[:p, :p]
            Pyx = Py[:, :p]
            covb += 1 / datapts2 * (np.kron(Px, Py) + comm_mat(pa, N * pa) @ np.kron(Pyx, Pyx.T))

        covb /= trials
        Wm = sp_pinv_with_retry(covb)
    else:
        # Use identity-based weighting
        Wm = np.eye(len(Eyy))

    ##################################################
    # Building the constant matrix for the LS problem
    ##################################################

    Ain = Aa - Aa @ L @ Ca
    OO = []
    temp = np.eye(n)
    for i in range(N):
        OO.append(Ca @ temp)
        temp = temp @ Ain
    OO = np.vstack(OO)

    # temporary variables
    M1 = np.zeros((na ** 2, ga ** 2))
    i = 0
    for j in range(ga):
        for k in range(ga):
            II = np.zeros((ga, ga))
            II[k, j] = 1
            # fun note: python-control's dlyaq does not work with non-symmetric Q's (Ga @ II @ Ga.T)... in Octave
            # this problem is solved by simply using the Sylvester equation with adapted inputs, even in the standard
            # dlyaq call. Here we have to use the solve_discrete_lyapunov function from scipy.linalg. The documentation
            # of octave-dlyaq does not mention this solution and it cost quite some time to find a function that does
            # the same.
            t1 = solve_discrete_lyapunov(Ain, Ga @ II @ Ga.T)
            M1[:, i] = t1.ravel(order="F")
            i += 1

    M2 = np.zeros((na ** 2, pa ** 2))
    i = 0
    for j in range(pa):
        for k in range(pa):
            II = np.zeros((pa, pa))
            II[k, j] = 1
            t2 = solve_discrete_lyapunov(Ain, Aa @ L @ II @ L.T @ Aa.T)
            M2[:, i] = t2.ravel(order="F")
            i += 1

    ##########################################
    # Single column ALS method
    ##########################################

    PSI = np.eye(pa)
    for i in range(N - 1):
        PSI = np.vstack([PSI, -Ca @ np.linalg.matrix_power(Ain, i) @ Aa @ L])

    OOtemp = np.kron(Ca, OO)
    PSItemp = np.kron(np.eye(pa), PSI)

    LH1 = OOtemp @ M1
    LH2 = OOtemp @ M2 + PSItemp

    if Rsym == 1:
        LHSsingc = np.hstack([LH1 @ symtran(ga), LH2 @ symtran(pa)])
    else:
        LHSsingc = np.hstack([LH1 @ symtran(ga), LH2[:, :pa ** 2:pa + 1]])

    M = np.kron(Ca, np.eye(n)) @ M1 @ symtran(ga)
    Mrank = np.linalg.matrix_rank(M, tol=0.0001)
    nr, nc = M.shape

    if nc > Mrank:
        print('WARNING: als_sdp_mrQ: The covariance estimates are not unique!\n'
              'Use trade-off curve to find lowest rank solution\n')

    if hasattr(estimator, 'Q') and hasattr(estimator, 'R'):
        Q0 = celldiag([estimator.Q, estimator.R])
    else:
        Q0 = np.eye(ga + pa)

    numax = 1
    Qest_cell = []
    Rest_cell = []

    lam_length = len(lam_vec) if lam_vec is not None else 1
    Phi = np.zeros(lam_length)
    Phi_tot = np.zeros(lam_length)
    Iter = np.zeros(lam_length)
    Iter_maxed = np.zeros(lam_length)
    trQ = np.zeros(lam_length)
    timespent = np.zeros(lam_length)

    if not sdp:
        # Solve least-squares problem without semi-definite constraints
        if lam_vec is not None:
            print('WARNING: als_sdp_mrQ: Cannot solve indefinite problem with penalty on trQ')

        time1 = process_time()
        X = np.linalg.inv(LHSsingc.T @ Wm @ LHSsingc) @ (LHSsingc.T @ Wm @ Eyy)
        timespent = process_time() - time1

        Qest_cell.append((symtran(ga) @ X[:ga * (ga + 1) // 2]).reshape((ga, ga), order="F"))
        if Rsym == 1:
            Rest_cell.append((symtran(pa) @ X[ga * (ga + 1) // 2:]).reshape((pa, pa), ordeR="F"))
        else:
            Rest_cell.append(np.diag(X[len(X) - pa:]))

        Phi = (Eyy - LHSsingc @ X).T @ Wm @ (Eyy - LHSsingc @ X)
        phi0 = Phi
        Phi_tot = Phi
        Iter = 0
        Iter_maxed = 0
        trQ = np.trace(Qest_cell[0] @ np.diag(trstates))
    else:
        # Solve ALS with semi-definite constraints
        # Scale weighting matrix and initial objective function
        res = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, 0, numax, pa, Rsym, trstates)
        QR0, phi0 = res[:2]
        Wm = Wm / phi0

        lam_iterator = lam_vec if lam_vec is not None else [0]

        # Loop through and solve for each value of lam_vec
        for i in range(len(lam_iterator)):
            lam = lam_iterator[i]

            time1 = process_time()
            res = sdp_QR_mrQ(LHSsingc, Eyy, Wm, Q0, lam, numax, pa, Rsym, trstates)
            timespent[i] = process_time() - time1

            QR1, phi, phi_tot, iter, nsq, iter_maxed = res[:6]
            Qest = QR1[:ga, :ga]
            Rest = QR1[ga:, ga:]
            Qest_cell.append((Qest + Qest.T) / 2)
            Rest_cell.append((Rest + Rest.T) / 2)
            Phi[i] = phi
            Phi_tot[i] = phi_tot
            Iter[i] = iter
            Iter_maxed[i] = iter_maxed
            trQ[i] = np.trace(Qest_cell[i] @ np.diag(trstates))

    # Form matrices of autocovariances from data and theory
    if Rsym == 1:
        Eyy_thry = LHSsingc @ np.hstack([symtranT(ga) @ Qest_cell[0].ravel(order="F"),
                                         symtranT(pa) @ Rest_cell[0].ravel(order="F")])
    else:
        Eyy_thry = LHSsingc @ np.hstack([symtranT(ga) @ Qest_cell[0].ravel(order="F"), np.diag(Rest_cell[0])])

    Eyy = Eyy.reshape(N * pa, pa, order="F")
    Eyy_thry = Eyy_thry.reshape(N * pa, pa, order="F")

    bhat_mtr = np.zeros((N, pa ** 2))
    bhat_mtr_thry = np.zeros((N, pa ** 2))
    for i in range(pa):
        for j in range(pa):
            for k in range(N):
                bhat_mtr[k, pa * i + j] = Eyy[pa * k + i, j]
                bhat_mtr_thry[k, pa * i + j] = Eyy_thry[pa * k + i, j]

    npts = datapts - start - N
    cov_bound = np.zeros((pa, pa))
    for i in range(pa):
        for j in range(pa):
            cov_bound[i, j] = np.sqrt((bhat_mtr[0, pa * i + i] * bhat_mtr[0, pa * j + j]) / npts)
    cov_bound = cov_bound.ravel(order="F")

    if plot_flag:
        n_cols = n_rows = pa
        fig = plt.figure()

        for i in range(pa ** 2):
            ax = fig.add_subplot(n_rows, n_cols, i + 1)

            ax.plot(range(N), bhat_mtr[:, i], '-*', range(N), bhat_mtr_thry[:, i], '-o', linewidth=2, label="Data")
            ax.plot([0, N], 2 * np.kron([[1, -1], [1, -1]], cov_bound[i]), 'k', label="Fit")
            if i == pa ** 2 - 1:
                ax.legend(bbox_to_anchor=(1.1, 1.05))
            if i <= pa:
                ax.set_title('Autocorrelation')
            if i > pa ** 2 - pa:
                ax.set_xlabel('Lag')

        plt.tight_layout()
        plt.show()

    return (Qest_cell, Rest_cell, trQ, Phi, Phi_tot, phi0, bhat_mtr, bhat_mtr_thry, cov_bound, Iter, Iter_maxed,
            timespent)
