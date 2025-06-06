# Generated with SMOP  0.41-beta
from ..libsmop import *


# als\inst\sdp_QR_mrQ.m

## SemiDefinite programming subroutine.  Uses Newton steps to satisfy
## the KKT conditions with log-barrier function


@function
def sdp_QR_mrQ(A=None, b=None, W=None, Q0=None, lam=None, nu=None, pa=None, Rsym=None, trstates=None, *args, **kwargs):
    varargin = sdp_QR_mrQ.varargin
    nargin = sdp_QR_mrQ.nargin

    # warning("off","Octave:singular-matrix-div")
    ntot = size(Q0, 1)
    # als\inst\sdp_QR_mrQ.m:7
    na = ntot - pa
    # als\inst\sdp_QR_mrQ.m:9
    ind = reshape(concat([arange(1, ntot ** 2)]), ntot, ntot)
    # als\inst\sdp_QR_mrQ.m:11
    ind = multiply((logical_and(ind, celldiag(cellarray([rand(na), rand(pa)])))), ind)
    # als\inst\sdp_QR_mrQ.m:12
    Az = zeros(size(A, 1), ntot ** 2)
    # als\inst\sdp_QR_mrQ.m:13
    # Reconstructs least squares matrix for R:
    if (pa != 0):
        LH2 = A(arange(), arange(dot(na, (na + 1)) / 2 + 1, end()))
        # als\inst\sdp_QR_mrQ.m:17
        if Rsym == 1:
            # For symmetric R matrix
            mtrx = symtranT(pa)
        # als\inst\sdp_QR_mrQ.m:20
        else:
            # For diagonal R matrix
            mtrx = zeros(pa, pa ** 2)
            # als\inst\sdp_QR_mrQ.m:23
            for i in arange(1, pa, 1).reshape(-1):
                mtrx[i, dot((i - 1), (pa + 1)) + 1] = 1
        # als\inst\sdp_QR_mrQ.m:25
        LH2 = dot(LH2, mtrx)
    # als\inst\sdp_QR_mrQ.m:28

    # Reconstructs least-squares matrix for Q:
    A = concat([dot(A(arange(), arange(1, dot(na, (na + 1)) / 2)), symtranT(na)), LH2])
    # als\inst\sdp_QR_mrQ.m:32
    idx = setdiff(ravel(ind), 0)
    # als\inst\sdp_QR_mrQ.m:33
    Az[arange(), idx] = A
    # als\inst\sdp_QR_mrQ.m:34
    alphmax = 1
    # als\inst\sdp_QR_mrQ.m:35

    ## Initialize
    tol = 1e-10
    # als\inst\sdp_QR_mrQ.m:38

    reltol = 1e-05
    # als\inst\sdp_QR_mrQ.m:39

    itertot = 0
    # als\inst\sdp_QR_mrQ.m:40
    maxiter = 100
    # als\inst\sdp_QR_mrQ.m:41

    Q = copy(Q0)
    # als\inst\sdp_QR_mrQ.m:42

    Qin = inv(Q)
    # als\inst\sdp_QR_mrQ.m:43
    con = 1
    # als\inst\sdp_QR_mrQ.m:44
    con1 = copy(con)
    # als\inst\sdp_QR_mrQ.m:45
    con2 = copy(con)
    # als\inst\sdp_QR_mrQ.m:46
    tol1 = copy(tol)
    # als\inst\sdp_QR_mrQ.m:47
    tol2 = copy(tol)
    # als\inst\sdp_QR_mrQ.m:48
    phi = 0
    # als\inst\sdp_QR_mrQ.m:49
    M = dot(dot(dot(2, Az.T), W), Az)
    # als\inst\sdp_QR_mrQ.m:50

    dtemp = dot(dot(dot(2, Az.T), W), b)
    # als\inst\sdp_QR_mrQ.m:51

    vlam_mtr = dot(lam, celldiag(cellarray([diag(trstates), zeros(pa)])))
    # als\inst\sdp_QR_mrQ.m:53
    vlam = ravel(vlam_mtr)
    # als\inst\sdp_QR_mrQ.m:54
    iter = 0
    # als\inst\sdp_QR_mrQ.m:55
    phi_old = 1000.0
    # als\inst\sdp_QR_mrQ.m:57
    nu_first = 1
    # als\inst\sdp_QR_mrQ.m:58
    while (nu > 1e-12):

        # Decrease penalty on log-barrier term (for sdp constraints)
        if iter > 0:
            if iter / maxiter < 0.75:
                nu = dot(nu, 0.75)
            # als\inst\sdp_QR_mrQ.m:63
            else:
                nu = dot(nu, (dot(0.95, iter))) / maxiter
        # als\inst\sdp_QR_mrQ.m:65
        iter = 1
        # als\inst\sdp_QR_mrQ.m:68
        phi = obj_tot(Q, Az, b, W, nu, vlam_mtr)
        # als\inst\sdp_QR_mrQ.m:69
        con1 = abs(phi - phi_old)
        # als\inst\sdp_QR_mrQ.m:70
        tol1 = tol + dot(reltol, abs(phi))
        # als\inst\sdp_QR_mrQ.m:71
        # Optimize for each penalty on log-barrier term
        while ((con1 > tol1 or con2 > tol2) and iter < maxiter):

            t1 = copy(cputime)
            # als\inst\sdp_QR_mrQ.m:74
            itertot = itertot + 1
            # als\inst\sdp_QR_mrQ.m:75
            iter = iter + 1
            # als\inst\sdp_QR_mrQ.m:76
            LH = M + dot(nu, kron(Qin, Qin))
            # als\inst\sdp_QR_mrQ.m:77
            RH = dot(M, ravel(Q)) - dtemp + vlam - dot(nu, ravel(Qin))
            # als\inst\sdp_QR_mrQ.m:78
            delQ = - reshape(dot(pinv(LH), RH), ntot, ntot)
            # als\inst\sdp_QR_mrQ.m:79
            # delQ = -reshape(LH\RH,ntot,ntot); # Can't calculate this way if Hessian is poorly conditioned
            ## Optimal Step Size
            alph = golden_section_Q_mrQ(0, alphmax, Q, delQ, Az, b, W, nu, vlam_mtr)
            # als\inst\sdp_QR_mrQ.m:82
            Qold = copy(Q)
            # als\inst\sdp_QR_mrQ.m:83
            Qnew = Qold + dot(dot(0.99, alph), delQ)
            # als\inst\sdp_QR_mrQ.m:84
            if (nu_first < 1 and min(eig(Qnew)) < 0):
                Qnew = copy(Qold)
            # als\inst\sdp_QR_mrQ.m:87
            # Update Q
            Q = copy(Qnew)
            # als\inst\sdp_QR_mrQ.m:90
            Q = (Q + Q.T) / 2
            # als\inst\sdp_QR_mrQ.m:91
            Qin = pinv(Q)
            # als\inst\sdp_QR_mrQ.m:92
            phi_old = copy(phi)
            # als\inst\sdp_QR_mrQ.m:94
            phi = obj_tot(Q, Az, b, W, nu, vlam_mtr)
            # als\inst\sdp_QR_mrQ.m:95
            con1 = abs(phi - phi_old)
            # als\inst\sdp_QR_mrQ.m:96
            con2 = norm(Q - Qold)
            # als\inst\sdp_QR_mrQ.m:97
            tol1 = tol + dot(reltol, abs(phi))
            # als\inst\sdp_QR_mrQ.m:98
            tol2 = tol + dot(reltol, norm(Q))
            # als\inst\sdp_QR_mrQ.m:99
            t2 = copy(cputime)
            # als\inst\sdp_QR_mrQ.m:100
            if t2 - t1 > 5:
                # Give warning if code is too slow
                warning('ALS iterations are really slow for some reason - you might want to check what is going on')

        # If solution with highest penalty on log-barrier term is infeasible, increase log-barrier penalty:
        if nu_first == 1:
            r, p = chol(Q, nargout=2)
            # als\inst\sdp_QR_mrQ.m:108
            if p == 0:
                nu_first = 0
            # als\inst\sdp_QR_mrQ.m:111
            else:
                # If solution is infeasible, increase nu
                nu = dot(nu, 100)
                # als\inst\sdp_QR_mrQ.m:114
                iter = 0
                # als\inst\sdp_QR_mrQ.m:115
                Q = copy(Q0)
                # als\inst\sdp_QR_mrQ.m:116
                Qin = inv(Q)
                # als\inst\sdp_QR_mrQ.m:117
                con = 1
                # als\inst\sdp_QR_mrQ.m:118
                con1 = copy(con)
                # als\inst\sdp_QR_mrQ.m:119
                con2 = copy(con)
                # als\inst\sdp_QR_mrQ.m:120
                tol1 = copy(tol)
                # als\inst\sdp_QR_mrQ.m:121
                tol2 = copy(tol)
                # als\inst\sdp_QR_mrQ.m:122
                phi = 0
            # als\inst\sdp_QR_mrQ.m:123
            # Make sure we can't get stuck in this loop
            if nu > 10000000000.0:
                error('sdp_QR_mrQ: Unable to find >0 solution')

    phi2 = obj_ls(Q, Az, b, W)
    # als\inst\sdp_QR_mrQ.m:132

    phi_uw = obj_ls(Q, Az, b, eye(size(W)))
    # als\inst\sdp_QR_mrQ.m:133

    nsq = dot(- nu, sum(log(eig(Q))))
    # als\inst\sdp_QR_mrQ.m:134

    iter_maxed = 0
    # als\inst\sdp_QR_mrQ.m:135
    if iter == maxiter:
        iter_maxed = 1
    # als\inst\sdp_QR_mrQ.m:137
